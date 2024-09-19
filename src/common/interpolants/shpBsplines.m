function shp = shpBsplines(n,x,~,xk,pMax)
%% SHPBSPLINES
% B-spline shape functions at local coordinate X for an edge with N nodes
% knot vector XK, spline order PMAX
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if nargin < 5                                                               % if order is not provided
    pMax=n-1;                                                               % choose p=n-1
end
if nargin < 4                                                               % if knot vector not provided
    xk=[-ones(1,pMax+1),ones(1,pMax+1)];                                    % create simple spline of order pMax
end

nStart=numel(xk)-1;                                                         % initial number of terms

tol=1e-9;                                                                   % tolerance for identifying identical knots

%% zeroth order
N=zeros(1,nStart);                                                          % shape functions
dN=zeros(1,nStart);                                                         % derivatives

ind=(x>=(xk(1:end-1)-tol))&(x<(xk(2:end)-tol));                             % Finding index of constant term in zeroth order spline
N(ind)=1;

for p=1:pMax                                                                % loop over orders
    i=1:numel(N)-1;
    
    D1=(xk(i+p)-xk(i));                                                     % demoninator, first term
    indZero=(abs(D1)<1e-10*max(abs(xk)));                                   % find zero values
    R1=( x-xk(i))./D1;
    R1(indZero)=0;                                                          % remove zero values
    
    D2=(xk(p+1+i)-xk(i+1));                                                 % denominator, second term
    indZero2=(abs(D2)<1e-10*max(abs(xk)));                                  % find zero values
    R2=( -x + xk(p+1+i) )./D2;
    R2(indZero2)=0;
    
    if p==pMax
        R3= p./D1;
        R3(indZero)=0;
        R4=-p./D2;
        R4(indZero2)=0;
        
        dN= R3.* N(i) + R4 .* N(i+1);
    end
    
    N= R1.* N(i) + R2 .* N(i+1);
    
end


shp=[N;dN];                                                                 % assemble shape functions and derivatives




