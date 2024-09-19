% kth derivative of jacobi polynomial
% https://en.wikipedia.org/wiki/Jacobi_polynomials#Derivatives

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
function jac = jacobidx(a,b,n,x,k)

% return first derivative if not specified
if nargin<5 
    k=1;
end
if (n-k)>=0
    jac=gamma(a+b+n+1+k)/2^k/gamma(a+b+n+1)*jacobi(a+k,b+k,n-k,x);
else
    jac=zeros(size(x));
end
