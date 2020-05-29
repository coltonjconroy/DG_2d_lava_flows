!*************************************************************************************
! Declaration of Precision of Variables
!*************************************************************************************

MODULE precisions

IMPLICIT NONE
SAVE

INTEGER, PARAMETER :: int_p  = SELECTED_INT_KIND(16)
INTEGER, PARAMETER :: real_p = SELECTED_REAL_KIND(16)

END MODULE precisions
!*************************************************************************************
!....Global variables
!*************************************************************************************
module global
use precisions
implicit none

integer(int_p) :: i, j, k, ii, p, nedges, nelems, nnodes, ndof, napts, nlpts
integer(int_p) :: nregions, i_count, j_count, k_count, irk, plot_n, nspts
integer(int_p) :: vid, wall_bc_type, limiter 
integer(int_p), allocatable :: ELEMnodes(:,:),ELEMedges(:,:)
integer(int_p), allocatable :: EDGEtype(:),EDGEelems(:,:),EDGEnodes(:,:)
integer(int_p), allocatable :: inlet_edges(:), outlet_edges(:)
integer(int_p), parameter :: debug = 0
integer(int_p), parameter :: outlet_bc_type = 0  !...0 = Df/Ds = 0; 1 = radiation bc
real(real_p), allocatable :: A(:,:,:), B(:,:,:), C(:,:), SC(:,:)
real(real_p), allocatable :: PHI(:,:), PHIedge(:,:,:)
real(real_p), allocatable :: PSI(:,:), PSIedge(:,:,:)
real(real_p), allocatable :: hu(:,:,:), hv(:,:,:), eta(:,:,:), mu(:,:)
real(real_p), allocatable :: rhs_hu(:,:,:), rhs_hv(:,:,:), rhs_eta(:,:,:)
real(real_p), allocatable :: ht(:,:,:), rhs_ht(:,:,:), dvdx(:,:,:), dudy(:,:,:)
real(real_p), allocatable :: ELEMarea(:), PHIedgef(:,:,:), EDGEpts(:,:,:)
real(real_p), allocatable :: ELEMxy(:,:), EDGEnormals(:,:), SRCEpts(:,:)
real(real_p), allocatable :: EDGElengths(:), Xnodes(:), Ynodes(:)
real(real_p), allocatable :: alpha(:,:), beta(:,:), ct(:), Z(:)
real(real_p), allocatable :: Zx(:), Zy(:), u_inlet(:), v_inlet(:)
real(real_p), allocatable :: Zcos(:), c_bry(:,:), PHInode(:,:),w(:,:,:)
real(real_p), allocatable :: dhdx(:,:),dhdy(:,:),dwdx(:,:),dwdy(:,:),b_thick(:,:)
real(real_p) :: t, dt, deltat, cfl, rho, T_int, hint, T_wall, T_crust, T_air
real(real_p) :: eff, b0, mu0, kt, vid_frame_rate, b_layer, alph
real(real_p) :: Amu, Bmu, Cmu, tau_yield
real(real_p), parameter :: pi = 3.1415926535897932d0, quart = 0.250d0
real(real_p), parameter :: g = 9.810d0, half = 0.50d0
real(real_p), parameter :: one_hf = 1.50d0
real(real_p), parameter :: mu_int = 0.050d0
real(real_p), parameter :: man_n  = 0.030d0
real(real_p), parameter :: sig = 5.670d-8          !...Boltzmann constant
real(real_p), parameter :: eHe = 7.50d+3           !...energy heat equivalent
real(real_p), parameter :: shc = 0.200d0           !...specific heat capacity of lava


end module global
!**************************************************************************
! Depth-integrated 2D lava flow model
! Performs time stepping and numerical integration
!
! input: A    - ndof x npts2d x 2   (Matrix) advection matrix (x-dir & y-dir)
!        B    - ndof x npts1d x 3   (Matrix) boundary matrix (3 edges)
!        C    - ndof x nptsL2       (Matrix) source matrix
!       PHI   - npts2d x ndof       (Matrix) triangular shape function
!   PHIedge   - npts1d x ndof x 3   (Matrix) edge shape function
!      PSI    - npts2d x 3          (Matrix) nodal element shape function
!   PSIedge   - npts1d x 3    x 3   (Matrix) nodal edge shape function
! ELEMnodes   - nelems x 3          (Matrix) element node connectivity
!  ELEMarea   - nelems x 1          (Vector) element areas
!    ELEMxy   - nelems x 4          (Matrix) xy Jacobians
! EDGEnormals - nedges x 2          (Matrix) edge normals
!  EDGEelems  - nedges x 2          (Matrix) edge element connectivity
!   EDGEtype  - nedges x 1          (Vector) edge type (1 == boundary)
! EDGElengths - nedges x 1          (Vector) edge length
! ELEMedges   - nelems x 3          (Matrix) elem edge connectivity
!         Z   - nnodes x 1          (Vector) thickness
!         t   -      1 x 1          (Scalar) local time
!         dt  -      1 x 1          (Scalar) advective time step
!
! Author: Colton J. Conroy
!         Columtia University
!         Lamont-Doherty Earth Observatory & Applied Physics and Applied Mathematics
!         August 2018 - May 2020
!
!**************************************************************************
program DI_lava_main

use precisions
use global
implicit  none

integer(int_p) :: RKorder, nrk, stages, new, NT, n, nframes
real(real_p) :: start, finish, Tend, time, check

!...user input (Tend, time-step)

write(*,*)'Welcome to the 2D depth-integrated lava flow model.'
write(*,*)'This program provides insight into unsteady depth-integrated lava flows.'
write(*,*)'Please insure all input files are ready to access...'
write(*,*)'Enter total time to run simulation and time step...'
read(*,*)Tend, dt

write(*,*)'Total simulation time', Tend, 'secs'
write(*,*)'Time step is', dt, 'secs'
!...pre-processes (read in data)

call read_dg_input()

call read_mesh_input()

call read_input_parameters()

!...explicit time steppers

RKorder = 4
stages  = 5
new = 1 !....New optimal time steppers? (1 = yes)
call rk_coefficients(RKorder,stages,new)
nrk = stages

!...initial condition

call read_initial_condition(nrk)

!....Video setup

plot_n  = 0
nframes = 0
if (vid == 1) then
    open(17,FILE='time.1d')
endif

!------------------------
! Main time stepping loop
!------------------------
call cpu_time(start)
NT = floor(Tend/dt)-1; time = 0.0d0;
do n = 0, NT
    do irk = 1, nrk
        t = time + ct(irk)*dt
        call rhs_dg_thermo(n)
        do i = 1, irk
            eta(:,:,irk+1) = eta(:,:,irk+1) + alpha(irk,i)*eta(:,:,i) &    ! thickness perturbation
                             + beta(irk,i)*dt*rhs_eta(:,:,i)
            hu(:,:,irk+1)  = hu(:,:,irk+1)  + alpha(irk,i)*hu(:,:,i)  &    ! hu
                             + beta(irk,i)*dt*rhs_hu(:,:,i)
            hv(:,:,irk+1)  = hv(:,:,irk+1)  + alpha(irk,i)*hv(:,:,i)  &    ! hv
                             + beta(irk,i)*dt*rhs_hv(:,:,i)
            ht(:,:,irk+1)  = ht(:,:,irk+1)  + alpha(irk,i)*ht(:,:,i)  &    ! ht
                             + beta(irk,i)*dt*rhs_ht(:,:,i)
        enddo
        if (limiter == 1) then
!            call slope_limiter() !<in progress>
        endif
    enddo
    hu(:,:,1)   = hu(:,:,nrk+1);    hu(:,:,2:nrk+1)  = 0.0d0; rhs_hu  = 0.0d0;
    hv(:,:,1)   = hv(:,:,nrk+1);    hv(:,:,2:nrk+1)  = 0.0d0; rhs_hv  = 0.0d0;
    ht(:,:,1)   = ht(:,:,nrk+1);    ht(:,:,2:nrk+1)  = 0.0d0; rhs_ht  = 0.0d0;
    eta(:,:,1)  = eta(:,:,nrk+1);   eta(:,:,2:nrk+1) = 0.0d0; rhs_eta = 0.0d0;
    w(:,:,1)    = w(:,:,nrk);       w(:,:,2:nrk)     = 0.0d0;
    dudy(:,:,1) = dudy(:,:,nrk);    dudy(:,:,2:nrk)  = 0.0d0;
    dvdx(:,:,1) = dvdx(:,:,nrk);    dvdx(:,:,2:nrk)  = 0.0d0;
    if (vid == 1) then
        call video(n,nframes,time)
    endif
    time  = time + dt
    check = maxval(hu(1,:,1))
    if (check > 10000.0d0) then
        write(*,*)'***** WARNING *****'
        write(*,*)'Numerical solution is blowing up.'
        write(*,*)'Execution will stop...'
        write(*,*)'Blow up occured at',time,'secs'
        goto 10
    endif
enddo
10 call cpu_time(finish)
!---------------
! write output
!--------------
call global_output()

write(*,*)'cpu time =',finish-start
write(*,*)'nframes =',nframes
stop
end program DI_lava_main
