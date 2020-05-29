!***********************************************************************
!
!     Subroutine read_initial_condition()
!
!     Reads in initial conditions and allocate main variables
!
!     Written by Colton J. Conroy (08-27-2018)
!
!***********************************************************************

subroutine read_initial_condition(nrk)

use precisions
use global
implicit none

integer(int_p) jj, kk, mesh, nrk, ndof_ic, nelems_ic, nnodes_ic
integer(int_p), parameter :: constant_h = 1

open(10,FILE='DI_lava_ic.txt')
read(10,'(A24)') mesh
read(10,*) ndof_ic, nelems_ic, nnodes_ic

!...check that ndof, nelems, and nnodes is consistent with other input files

if (ndof_ic == ndof) then
    if (nelems_ic == nelems) then
        if (nnodes_ic == nnodes) then

            allocate(   hu(ndof,nelems,nrk+1) )
            allocate(   hv(ndof,nelems,nrk+1) )
            allocate(rhs_hu(ndof,nelems,nrk)  )
            allocate(rhs_hv(ndof,nelems,nrk)  )
            allocate(     Z(nnodes)           )
            allocate(    Zx(nnodes)           )
            allocate(    Zy(nnodes)           )
            allocate(   Zcos(nnodes)          )
            allocate(    u_inlet(nnodes)      )
            allocate(    v_inlet(nnodes)      )
            allocate( eta(ndof,nelems,nrk+1)  )
            allocate(  ht(ndof,nelems,nrk+1)  )
            allocate(  mu(ndof,nelems)        )
            allocate(rhs_eta(ndof,nelems,nrk) )
            allocate( rhs_ht(ndof,nelems,nrk) )
            allocate( dvdx(ndof,nelems,nrk)   )
            allocate( dudy(ndof,nelems,nrk)   )
            allocate( dhdx(ndof,nelems)       )
            allocate( dhdy(ndof,nelems)       )
            allocate( dwdx(ndof,nelems)       )
            allocate( dwdy(ndof,nelems)       )
            allocate(    w(ndof,nelems,nrk)   )
            allocate( b_thick(ndof,nelems)    )

            hu(:,:,:)      = 0.0d0;     hv(:,:,:)  = 0.0d0;
            rhs_hu(:,:,:)  = 0.0d0; rhs_hv(:,:,:)  = 0.0d0;
            z(:)           = 0.0d0;         zx(:)  = 0.0d0;
            zy(:)          = 0.0d0;    eta(:,:,:)  = 0.0d0;
            mu(:,:)        = 0.0d0;     ht(:,:,:)  = 0.0d0;
            rhs_eta(:,:,:) = 0.0d0; rhs_ht(:,:,:)  = 0.0d0;
            dvdx(:,:,:)    = 0.0d0;   dudy(:,:,:)  = 0.0d0;
            Zcos(:)        = 0.0d0;   dhdx(:,:)    = 0.0d0;
            dhdy(:,:)      = 0.0d0;   b_thick(:,:) = 0.0d0;
            dwdx(:,:)      = 0.0d0;   dwdy(:,:)    = 0.0d0;
            w(:,:,:)       = 0.0d0;

            do i = 1, nelems
                do j = 1, ndof
                    read(10,*) ii,hu(j,i,1),hv(j,i,1)
                    eta(j,i,1) = 0.0d0
                    if (j == 1) then
                        mu(j,i)   = mu_int
                        ht(j,i,1) = hint*T_int
                    else
                        mu(j,i)   = 0.0d0
                        ht(j,i,1) = 0.0d0
                    endif
                enddo
            enddo

            do i = 1, nnodes
                read(10,*) ii, Z(i), Zx(i), Zy(i), Zcos(i)
            enddo
            if (constant_h == 1) then
                Z(:) = hint
            endif
            do i = 1, nnodes
                read(10,*) ii, u_inlet(i), v_inlet(i)
            enddo

        else
            write(*,*)'input files do not match'
            write(*,*)'  execution will stop.  '
            write(*,*)'ndof=',ndof
            write(*,*)'nelems=',nelems
            write(*,*)'nnodes=',nnodes
            write(*,*)'ndof=',ndof_ic
            write(*,*)'nelems=',nelems_ic
            write(*,*)'nnodes=',nnodes_ic
            stop
        endif
    endif
endif

close(10)
end subroutine read_initial_condition
