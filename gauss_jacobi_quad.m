function [X,W] = gauss_jacobi_quad(n,a,b)

%--------------------------------------------------------------------------
%
%  [X,W] = gauss_jacobi_quad(n,a,b)
%   
%  Returns the n Gauss-Jacobi quadrature points, X (in ascending order) and 
%  their respective weights W to approximate an integral of the form
%
%  int( (1-x)^a*(1+x)^b*f(x),x,-1,1 ). 
%
%--------------------------------------------------------------------------
%
%  Input:
%  ------
%       n: number of Gauss points (required).
%    a, b: exponents on weight function (1-x)^a and (1+x)^b (optional).
%          Default is a=b=0, which corresponds to standard Gauss-Legendre
%          quadrature.
%
%  Output:
%  -------
%    X:  1 by n row vector of the n Gauss-Jacobi points
%    W:  1 by n row vector of the respective weights
%
%--------------------------------------------------------------------------
%  NOTES: (1) An n-point Gauss-Jacobi quadrature will integrate all
%             polynomials f(x) up to degree 2*n-1 exactly.
%         (2) The special case a = b = 0 corresponds to Gauss-Legendre
%             quadrature.
%--------------------------------------------------------------------------
%
%  Written by Ethan Kubatko, July 4, 2011
%
%--------------------------------------------------------------------------

if nargin == 1
    a = 0; b = 0;
end

% Construct the N-th Jacobi polynomial
%-------------------------------------

P    = jacobi_poly(n,a,b);
dPdx = polyder(P);

% Obtain the n Gauss-Jacobi points and sort in ascending order
%-------------------------------------------------------------

X = sort(roots(P));

% Compute the corresponding n weights
%------------------------------------

W = gamma(a+n+1)*gamma(b+n+1)*2^(2*n+a+b+1)*factorial(n)./(...
    gamma(n+a+b+1)*(1-X.^2).*(polyval(dPdx,X)*2^n*factorial(n)/(-1)^n).^2);
end