!--------------------------------------------------------------------
!
!                   Read in model parmeters
!                 Written by Colton J. Conroy
!                       @ 480 Convent
!                         11.26.18
!
!-------------------------------------------------------------------
subroutine read_input_parameters()

use precisions
use global
implicit none

integer(int_p) :: title

!..open file

open(10,FILE='parameter_input.txt')
read(10,'(A24)') title

!..boundary conditions

read(10,*)wall_bc_type        ! 0 = no slip, 1 = no normal flow

!..dynamic parameters

read(10,*)rho
read(10,*)hint
read(10,*)b_layer
read(10,*)alph
read(10,*)Amu
read(10,*)Bmu
read(10,*)Cmu
read(10,*)tau_yield

!..thermal parameters

read(10,*)T_int
read(10,*)T_wall
read(10,*)T_crust
read(10,*)T_air
read(10,*)eff                  !...effusivity constant
read(10,*)b0                   !...experimental reference temperature
read(10,*)mu0                  !...viscosity @ infinite temperature
read(10,*)kt                   !...thermal conductivity

!..miscellaneous

read(10,*)vid
read(10,*)vid_frame_rate
read(10,*)limiter

close(10)

return
end subroutine read_input_parameters
