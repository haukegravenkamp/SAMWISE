function Pkl=proriol(k,l,x1,x2)
% Proriol polynomial of order k,l on the standard square
% Derived from the Jacobi polynomials

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
Pkl = jacobi(0,0,k,x1).*((1-x2)/2).^k.*jacobi(2*k+1,0,l,x2);


