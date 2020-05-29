!***********************************************************************
!
!     Subroutine read_dg_input()
!
!     Reads in dg matricies and basis functions
!
!     Written by Colton J. Conroy (08-24-2018)
!
!***********************************************************************

subroutine read_dg_input()

use precisions
use global
implicit none

integer(int_p) dg

open(10,FILE='dg_input.txt')
read(10,'(A24)') dg
read(10,*) ndof, napts, nlpts, nspts, p

allocate( A(ndof,napts,2)       )
allocate( B(ndof,nlpts,3)       )
allocate( C(ndof,napts)         )
allocate( SC(ndof,nspts)        )
allocate( PHI(napts,ndof)       )
allocate( PHInode(3,ndof)       )
allocate( PHIedge(nlpts,ndof,3) )
allocate( PSI(napts,3)          )
allocate( PSIedge(nlpts,3,3)    )
allocate(PHIedgef(nlpts,ndof,3) )

A(:,:,:) = 0.0d0; B(:,:,:)       = 0.0d0;  C(:,:) = 0.0d0;
PHI(:,:) = 0.0d0; PHIedge(:,:,:) = 0.0d0; SC(:,:) = 0.0d0;
PSI(:,:) = 0.0d0; PSIedge(:,:,:) = 0.0d0; PHInode(:,:) = 0.0d0;
PHIedgef(:,:,:) = 0.0d0

do i = 1, 2
    do j = 1, napts
        do k = 1, ndof
            read(10,*) A(k,j,i)
        enddo
    enddo
enddo

do i = 1, 3
    do j = 1, nlpts
        do k = 1, ndof
            read(10,*) B(k,j,i)
        enddo
    enddo
enddo

do i = 1, napts
    do j = 1, ndof
        read(10,*)C(j,i)
    enddo
enddo

do i = 1, nspts
    do j = 1, ndof
        read(10,*)SC(j,i)
    enddo
enddo

do i = 1, ndof
    do j = 1, napts
        read(10,*)PHI(j,i)
    enddo
enddo

do i = 1, ndof
    do j = 1, 3
        read(10,*)PHInode(j,i)
    enddo
enddo

do i = 1, 3
    do j = 1, ndof
        do k = 1, nlpts
            read(10,*) PHIedge(k,j,i)
        enddo
    enddo
enddo

do i = 1, 3
    do j = 1, napts
        read(10,*)PSI(j,i)
    enddo
enddo

do i = 1, 3
    do j = 1, 3
        do k = 1, nlpts
            read(10,*) PSIedge(k,j,i)
        enddo
    enddo
enddo

close(10)

do i = 1, 3
    do j = 1, ndof
        do k = 1, nlpts
            PHIedgef(k,j,i) = PHIedge(nlpts+1-k,j,i)
        enddo
    enddo
enddo

end subroutine read_dg_input
