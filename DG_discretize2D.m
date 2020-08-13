function [M,A,B,P,PHI,C,PSI,na_Gauss,nl_Gauss,S,Pn,EDGEpts,SRCEpts,Ps,Ws,ELEMpts] = DG_discretize2D(p)
%-----------------------------------------------------------------------------------
%   function [M,A,B,P,PHI,S,PSI,na_Gauss,nl_Gauss] = DG_discretize2D(p)
%-----------------------------------------------------------------------------------
%   DG master element matrices
%
%       p -  order of element for 2D DG spatial discretization
%       M -  mass matrix of size ndof x ndof
%       A - advection matrix of size ndof x na_Gauss x 2.
%             where m is the number of degrees of freedom for nonlinear
%             elements, and n is the fewest number of Gauss points needed to
%             integrate an area integral exactly and the first and second
%             pages of the matrix refer to the master element coordinate
%             system
%       B -  matrix of size m x n x 3.
%             holds the values of the edge integrals for triangular elements
%             where pages 1, 2, 3 correspond to edge integrals for local edges
%             1, 2 and 3, respectively
%       P -  matrix of size L x n
%             where L is length of the vector containing area Gauss points
%             n is number of polynomials for the so called Dubiner basis
%             Holds the values of the Dubiner basis functions evaluated at
%             the area Gauss points.
%       PHI - PHI(#).edge where # is 1,2,3.
%               holds the basis functions over local edges 1, 2, 3
%               respectively, based on edge Gauss values.
%       S -  source function matrix of size ndof x na_Gauss
%
%   Written by Angela Nappi. 
%   Updated by Colton J. Conroy
%-----------------------------------------------------------------------------------
ndof = (p+1)*(p+2)/2; 

% Cubature rule for area integral must integrate a polynomial of degree 2*p exactly
[Xa,Ya,Wa] = triangle_quad(2*p); close  % area integrals
[Xs,Ys,Ws] = triangle_quad(2*(p+2)); close  % area integrals
% Dubiner triangle basis functions evaluated at area Gauss points.
[P,dPdx,dPdy] = triangle_basis(p,Xa,Ya); 
Ps = triangle_basis(p,Xs,Ys); 
% Basis evaluated at nodes 
Xn = zeros(3,1); Yn = zeros(3,1);
Xn(1,1) = -1; Yn(1,1) = -1;       % node 1 
Xn(2,1) =  1; Yn(2,1) = -1;       % node 2 
Xn(3,1) = -1; Yn(3,1) =  1;       % node 3 
Pn = triangle_basis(p,Xn,Yn); 
% 1D Gauss quadrature rule for line integrals using a Legendre basis
[xl,wl] = gauss_jacobi_quad(p+1,0,0); 
na_Gauss = length(Xa);
nl_Gauss = length(xl);
ns_Gauss = length(Xs);
% calculate mass matrix
M = zeros(ndof,ndof);
for i = 1 : ndof
    for j = 1 : ndof
        for k = 1 : na_Gauss
        M(i,j) = M(i,j) + Wa(k) * P(k,i) * P(k,j);
        end
    end
end

% Calculate advection matrix & source matrix
A =zeros(ndof,na_Gauss,2);
C = zeros(ndof,na_Gauss);
for i = 1 : ndof
    for j = 1 : na_Gauss % no. of area Gauss points
        A(i,j,1) = A(i,j,1) + Wa(j) * dPdx(j,i);
        A(i,j,2) = A(i,j,2) + Wa(j) * dPdy(j,i);
        C(i,j)   = C(i,j)   + Wa(j) * P(j,i);
    end
end

S = zeros(ndof,ns_Gauss);
for i = 1 : ndof
    for j = 1 : ns_Gauss % no. of area Gauss points
        S(i,j) = S(i,j)   + Ws(j) * Ps(j,i);
    end
end
% Calculate edge integral matrices
B = zeros(ndof,nl_Gauss,3); 
E1xi1 = zeros(nl_Gauss,1); E1xi2 = zeros(nl_Gauss,1);
E2xi1 = zeros(nl_Gauss,1); E2xi2 = zeros(nl_Gauss,1);
E3xi1 = zeros(nl_Gauss,1); E3xi2 = zeros(nl_Gauss,1);
for j = 1 : nl_Gauss 
    E1xi1(j) = -xl(j);  
    E1xi2(j) = xl(j);  
    E2xi1(j) = -1;  
    E2xi2(j) = -xl(j); 
    E3xi1(j) = xl(j);  
    E3xi2(j) = -1; 
end

EDGEpts = zeros(nl_Gauss,2,3);
EDGEpts(:,1,1) = E1xi1;
EDGEpts(:,2,1) = E1xi2; 
EDGEpts(:,1,2) = E2xi1;
EDGEpts(:,2,2) = E2xi2;
EDGEpts(:,1,3) = E3xi1;
EDGEpts(:,2,3) = E3xi2;

ELEMpts        = zeros(na_Gauss,2);
ELEMpts(:,1)   = Xa;
ELEMpts(:,2)   = Ya;

SRCEpts        = zeros(ns_Gauss,2);
SRCEpts(:,1)   = Xs;
SRCEpts(:,2)   = Ys; 

[PHI(1).edge] = triangle_basis(p,E1xi1,E1xi2); % triangle basis function evaluated at Gauss points for local edge 1
    P1 = PHI(1).edge';
[PHI(2).edge] = triangle_basis(p,E2xi1,E2xi2); % "  " for local edge 2
    P2 = PHI(2).edge';
[PHI(3).edge] = triangle_basis(p,E3xi1,E3xi2); % "  " for local edge 3
    P3 = PHI(3).edge';
for i = 1 : ndof
   for j = 1 : nl_Gauss
       B(i,j,1) = B(i,j,1) + wl(j) * P1(i,j);
       B(i,j,2) = B(i,j,2) + wl(j) * P2(i,j);
       B(i,j,3) = B(i,j,3) + wl(j) * P3(i,j);
   end
end

% Wind interpolation structure PSI
for i = 1 : na_Gauss
    PSI(i).elem(1) = -0.5*(Xa(i)+Ya(i));
    PSI(i).elem(2) = 0.5*(1 + Xa(i));
    PSI(i).elem(3) = 0.5*(1 + Ya(i));
end

for j = 1 : nl_Gauss
    PSI(j).edge1(1) = -0.5*(E1xi1(j)+E1xi2(j));
    PSI(j).edge1(2) = 0.5*(1 + E1xi1(j));
    PSI(j).edge1(3) = 0.5*(1 + E1xi2(j));
    
    PSI(j).edge2(1) = -0.5*(E2xi1(j)+E2xi2(j));
    PSI(j).edge2(2) = 0.5*(1 + E2xi1(j));
    PSI(j).edge2(3) = 0.5*(1 + E2xi2(j));
    
    PSI(j).edge3(1) = -0.5*(E3xi1(j)+E3xi2(j));
    PSI(j).edge3(2) = 0.5*(1 + E3xi1(j));
    PSI(j).edge3(3) = 0.5*(1 + E3xi2(j));
end

% PSI(1).edge1 = -0.5.*(E1xi1+E1xi2);
% PSI(2).edge1 = 0.5.*(ones(nl_Gauss,1) + E1xi1);
% PSI(3).edge1 = 0.5.*(ones(nl_Gauss,1) + E1xi2);

% PSI(1).edge2 = -0.5.*(E2xi1+E2xi2);
% PSI(2).edge2 = 0.5.*(ones(nl_Gauss,1) + E2xi1);
% PSI(3).edge2 = 0.5.*(ones(nl_Gauss,1) + E2xi2);

% PSI(1).edge3 = -0.5.*(E3xi1+E3xi2);
% PSI(2).edge3 = 0.5.*(ones(nl_Gauss,1) + E3xi1);
% PSI(3).edge3 = 0.5.*(ones(nl_Gauss,1) + E3xi2);

end % function













    
    