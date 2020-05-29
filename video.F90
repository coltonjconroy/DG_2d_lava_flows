!--------------------------------------------------------------------
!
!                          Output for Video
!                     Written by Colton J. Conroy
!                           @ the C.H.I.L
!                               3.21.13
!
!---------------------------------------------------------------------
subroutine video(n,nframes,time)

use precisions
use global
implicit none

integer(int_p) :: n, REGION, nframes
real(real_p) :: time
character(50) ::  xstring, estring, ystring, tstring

if (plot_n == n) then
    nframes = nframes + 1
    write(estring,'(i10)')plot_n+1
    open(10,FILE=estring)
    do i = 1, nelems
        write(10,*)eta(:,i,1)
    enddo
    close(10)
    write(xstring,'(i10)')plot_n
    write(ystring,'(i10)')plot_n+2
    !....x-component of momentum
    open(11,FILE=xstring)
    do i = 1, nelems
        write(11,*)hu(:,i,1)
    enddo
    close(11)
    !....y-component of momentum
    open(12,FILE=ystring)
    do i = 1, nelems
        write(12,*)hv(:,i,1)
    enddo
    close(12)
    write(tstring,'(i10)')plot_n+3
    open(13,FILE=tstring)
    do i = 1, nelems
        write(13,*)ht(:,i,1)
    enddo
    close(13)
    write(17,*)time
    plot_n = plot_n + floor(vid_frame_rate/dt)
endif

return
end SUBROUTINE video
