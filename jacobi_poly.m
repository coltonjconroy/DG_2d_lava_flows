function [p] = jacobi_poly(n,a,b)

%--------------------------------------------------------------------------
%
%  [p] = jacobi_poly(n,a,b)
%   
%  Returns the n-th order Jacobi polynomial of weight (a,b), P_n^(a,b), 
%  represented as a row vector containing the coefficients ordered by 
%  descending powers. For example, the polynomial
% 
%  P_3^(0,0) = 5/2*x^3 - 3/2*x 
%
%  would be represented as
%
%  p = [5/2, 0, -3/2, 0 ].
%
%  This type of polynomial representation can be used with the various Mat-
%  lab polynomial functions, e.g., polyval, polyder, polyint, etc. 
%
%--------------------------------------------------------------------------
%
%  Input:
%  ------
%       n: degree of Jacobi polynomial.
%    a, b: exponents on weight function (1-x)^a and (1+x)^b
%
%  Output:
%  -------
%    p:  1 by n+1 row vector of the Jacobi polynomial coefficients ordered
%        by descending powers.
%
%--------------------------------------------------------------------------
%  NOTE: (1) The special case a = b = 0 corresponds to the n-th degree 
%            Legendre polynomial.
%--------------------------------------------------------------------------
%
%  Written by Ethan Kubatko, October 7, 2011
%
%--------------------------------------------------------------------------

p = zeros([1,n+1]);
for k = 0:n
    poly = 1;
    for i = 1:k
        poly = conv(poly,[1/2,-1/2]);
    end
    p(n+1-k:n+1) = p(n+1-k:n+1) + factorial(n)/(factorial(k)*...
                        factorial(n-k))*gamma(a+b+n+k+1)/gamma(a+k+1)*poly;
end
p = gamma(a+n+1)/(gamma(a+b+n+1)*factorial(n))*p;
end