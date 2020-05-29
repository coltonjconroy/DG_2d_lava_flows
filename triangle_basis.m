function [P,dPdx,dPdy] = triangle_basis(p,X,Y)

%--------------------------------------------------------------------------
%
%   P = triangle_basis(p,X,Y)
%
%   Evaluates the first n = (p+1)(p+2)/2 polynomials for the so called 
%   Dubiner basis for the triangle at the points X, Y. 
%
%--------------------------------------------------------------------------
%
%  Input:
%  ------
%    p:  highest degree of polynomial to evaluate     
%    X:  vector containing the x values to be evaluated
%    Y:  vector containing the y values to be evaluated
%
%  Output:
%  -------
%
%    The matrix P is ordered as follows:
%
%         P0(X(1),Y(1))    P1(X(1),Y(1))  .  .  Pn(X(1),Y(1))
%         P0(X(2),Y(2))    P1(X(2),Y(2))  .  .  Pn(X(2),Y(2))
%              .                .        .           .
%              .                .           .        .
%         P0(X(L),Y(L))    P1(X(L),Y(L))  .  .  Pn(X(L),Y(L))
%
%    Note: size(P) = L by (p+1)(p+2)/2 where L = length(X)
%
%--------------------------------------------------------------------------
%
%  Written by Ethan Kubatko
%
%--------------------------------------------------------------------------

% Construct the basis functions and their derivatives

syms x y z

for i=0:p    
    Phi_a(i+1,1) = expand(1/(2^i*factorial(i))*diff((z^2-1)^i,i,z));    
end

z = 2*(1+x)/(1-y) - 1;
Phi_a = eval(Phi_a);    

m = 1;
for j=0:p
    for i=0:p-j
        jacobi(i+1,j+1) = 1/((-1)^j*2^j*factorial(j)*(1-y)^(2*i+1)*(1+y)^0)*...
        diff((1-y)^(2*i+1)*(1+y)^0*(1-y^2)^j,j,y);
        Phi_b(i+1,j+1) = ((1-y)/2)^i*jacobi(i+1,j+1);
        Phi(i+1,j+1) = simplify(expand(Phi_a(i+1)*Phi_b(i+1,j+1)));
        m = m + 1;
    end
end

m = 1;
for j=0:p
    for i=0:j
        jj = j - i;
        PHI(m,1) = Phi(i+1,jj+1);
        m = m + 1;
    end
end
PHI;
dPHIdx = diff(PHI,x);
dPHIdy = diff(PHI,y);

% Evaluate the basis functions at the given points

for i=1:length(X)    
    x = X(i);
    y = Y(i);    
    P(i,:) = eval(PHI.');
    dPdx(i,:) = eval(dPHIdx.');
    dPdy(i,:) = eval(dPHIdy.');
end










