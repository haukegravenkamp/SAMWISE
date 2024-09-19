function Pkl=prorioldx2(k,l,x1,x2)
%% Derivative of Proriol polynomial of order k,l on the standard square

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
Pkl=jacobi(0,0,k,x1).*(...                                   % unchanged
    ((-k/2)*((1-x2)/2).^(k-1).*jacobi(2*k+1,0,l,x2)...       % product rule for x2
    + ((1-x2)/2).^k.*jacobidx(2*k+1,0,l,x2,1)) );


