!-----------------------------------------------------------------------------------
!
!                 Output Subroutine
!              Written by Colton J. Conroy
!                      @ LDEO
!                      7.13.18
!
!------------------------------------------------------------------------------------
subroutine global_output()

use precisions
use global
implicit none

!....Write output for visualization program

open(10,file='output.txt')

write(10,*)nelems
write(10,*)ndof
write(10,*)nnodes

do i = 1, nelems
    do j = 1, ndof
        write(10,*)hint + eta(j,i,1)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)hu(j,i,1)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)hv(j,i,1)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)ht(j,i,1)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)mu(j,i)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)dvdx(j,i,1)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)dudy(j,i,1)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)-dhdx(j,i)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)-dhdy(j,i)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)b_thick(j,i)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)w(j,i,1)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)dwdx(j,i)
    enddo
enddo

do i = 1, nelems
    do j = 1, ndof
        write(10,*)dwdy(j,i)
    enddo
enddo

close(10)

return
end subroutine global_output

