[M,A,B,P,PHI,C,PSI,na,nl,Sc,Pn,EDGEpts,SRCEpts,Ps,Ws] = DG_discretize2D(p);
A(:,:,1) = M\A(:,:,1);
A(:,:,2) = M\A(:,:,2);
B(:,:,1) = M\B(:,:,1);
B(:,:,2) = M\B(:,:,2);
B(:,:,3) = M\B(:,:,3);
C  = M\C; 
Sc = M\Sc;
PHIedge = zeros(nl,ndof,3);
PHIedge(:,:,1) = PHI(1).edge;
PHIedge(:,:,2) = PHI(2).edge;
PHIedge(:,:,3) = PHI(3).edge;
sizeM = zeros(ndof,nl);
Psi = zeros(na,3);
for i = 1:na
    Psi(i,:) = PSI(i).elem;
end
PSIedge = zeros(nl,3,3);
for i = 1:nl
    PSIedge(i,:,1) = PSI(i).edge1;
    PSIedge(i,:,2) = PSI(i).edge2;
    PSIedge(i,:,3) = PSI(i).edge3;
end
EDGElengths = zeros(nedges,1); EDGEtype = zeros(nedges,1);
EDGEelems   = zeros(nedges,2); EDGEnormals = zeros(nedges,2);
EDGEnodes   = zeros(nedges,2);
ELEMnodes   = zeros(nelems,3); ELEMarea = zeros(nelems,1);
ELEMxy      = zeros(nelems,4); ELEMedges = zeros(nelems,3);
for i = 1:nedges
    EDGElengths(i)   = EDGES(i).length;
    EDGEtype(i)      = EDGES(i).type;
    if EDGEtype(i) == 1
        EDGEelems(i,1)   = EDGES(i).elems(1);
        EDGEelems(i,2)   = 0;
    else
        EDGEelems(i,1)   = EDGES(i).elems(1);
        EDGEelems(i,2)   = EDGES(i).elems(2);
    end
    EDGEnormals(i,1) = EDGES(i).normal(1);
    EDGEnormals(i,2) = EDGES(i).normal(2);
    EDGEnodes(i,1)   = EDGES(i).nodes(1);
    EDGEnodes(i,2)   = EDGES(i).nodes(2);
end
for i = 1:nelems;
    ELEMarea(i)    = ELEMS(i).area;
    ELEMxy(i,1)    = ELEMS(i).y12;
    ELEMxy(i,2)    = ELEMS(i).x21;
    ELEMxy(i,3)    = ELEMS(i).x13;
    ELEMxy(i,4)    = ELEMS(i).y31;
    ELEMedges(i,1) = ELEMS(i).edges(1);
    ELEMedges(i,2) = ELEMS(i).edges(2);
    ELEMedges(i,3) = ELEMS(i).edges(3);
    ELEMnodes(i,1) = ELEMS(i).nodes(1);
    ELEMnodes(i,2) = ELEMS(i).nodes(2);
    ELEMnodes(i,3) = ELEMS(i).nodes(3);
end

%clear ELEMS EDGES PSI PHI