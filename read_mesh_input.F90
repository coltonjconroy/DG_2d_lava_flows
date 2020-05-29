!***********************************************************************
!
!     Subroutine read_mesh_input()
!
!     Reads in simplex mesh information
!
!     Written by Colton J. Conroy (08-23-2018)
!
!***********************************************************************

subroutine read_mesh_input()

use precisions
use global
implicit none

integer(int_p) mesh

open(10,FILE='mesh_input.txt')
read(10,'(A24)') mesh
read(10,*) nelems, nedges, nnodes

allocate( ELEMnodes(nelems,3)   )
allocate( ELEMarea(nelems)      )
allocate( ELEMxy(nelems,4)      )
allocate( ELEMedges(nelems,3)   )
allocate( EDGEnormals(nedges,2) )
allocate( EDGEelems(nedges,2)   )
allocate( EDGEtype(nedges)      )
allocate( EDGElengths(nedges)   )
allocate( EDGEnodes(nedges,2)   )
allocate( Xnodes(nnodes)        )
allocate( Ynodes(nnodes)        )
allocate( inlet_edges(nedges)   )
allocate( outlet_edges(nedges)  )
allocate( c_bry(nelems,2)       )
allocate( SRCEpts(nspts,2)      )
allocate( EDGEpts(nlpts,2,3)      )

ELEMnodes(:,:)   = 0;     ELEMarea(:)    = 0.0d0;
ELEMxy(:,:)      = 0.0d0; ELEMedges(:,:) = 0;
EDGEnormals(:,:) = 0.0d0; EDGEelems(:,:) = 0;
EDGEtype(:)      = 0;     EDGElengths(:) = 0.0d0;
EDGEnodes(:,:)   = 0.0d0; Xnodes(:)      = 0.0d0;
Ynodes(:)        = 0.0d0; inlet_edges(:) = 0.0d0;
outlet_edges(:)  = 0.0d0; c_bry(:,:)     = 0.0d0;
SRCEpts(:,:)     = 0.0d0; EDGEpts(:,:,:) = 0.0d0;

do i = 1, nelems
    read(10,*) ii,ELEMnodes(i,1),ELEMnodes(i,2),ELEMnodes(i,3)
enddo

do i = 1, nelems
    read(10,*) ii,ELEMarea(i)
enddo

do i = 1, nelems
    read(10,*) ii,ELEMxy(i,1),ELEMxy(i,2),ELEMxy(i,3),ELEMxy(i,4)
enddo

do i = 1, nelems
    read(10,*) ii,ELEMedges(i,1),ELEMedges(i,2),ELEMedges(i,3)
enddo

do i = 1, nedges
    read(10,*) ii,EDGEnormals(i,1),EDGEnormals(i,2)
enddo

do i = 1, nedges
    read(10,*) ii,EDGEelems(i,1),EDGEelems(i,2)
enddo

do i = 1, nedges
    read(10,*) ii,EDGEtype(i)
enddo

do i = 1, nedges
    read(10,*) ii,EDGElengths(i)
enddo

do i = 1,nedges
    read(10,*) ii, EDGEnodes(i,1), EDGEnodes(i,2)
enddo

do i = 1, nnodes
    read(10,*)ii,Xnodes(i),Ynodes(i)
enddo

do i = 1, nedges
    read(10,*)ii, inlet_edges(i)
enddo

do i = 1, nedges
    read(10,*)ii, outlet_edges(i)
enddo

do i = 1, nelems
    read(10,*)ii, c_bry(i,1), c_bry(i,2)
enddo

do i = 1, nspts
    read(10,*)ii,SRCEpts(i,1),SRCEpts(i,2)
enddo

do j = 1, 3
    do i = 1, nlpts
        read(10,*)ii,EDGEpts(i,1,j),EDGEpts(i,2,j)
    enddo
enddo

close(10)

return
end subroutine read_mesh_input
