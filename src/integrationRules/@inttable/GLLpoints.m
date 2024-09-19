function GaussLobatto = GLLpoints(N)
N = N-1;

% Computes the Legendre-Gauss-Lobatto nodes, weights 

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N)/N)';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
             
end

w=2./(N*N1*P(:,N1).^2);

x = -x;
d = zeros(1,N1);
for jj = 1:N1
    pts = x;
    b = x(jj) - pts;
    b(jj) = 1;
    d(jj) = 1/prod(b);
end

GaussLobatto = struct('xi', x', 'wgt', w', 'prd',d);
