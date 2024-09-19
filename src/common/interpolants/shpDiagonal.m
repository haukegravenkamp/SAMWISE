function shp = shpDiagonal(n,eta,inttable)
%% SHPDIAGONAL
% "Diagonal" shape functions
% Shape functions used in the "Diagonal SBFEM"
% These functions are polynomials of order 2n-1 (n: number of nodes)
% The polynonials have Kronecker-Delta property at the nodes.
% The derivatives vanish at the nodes.
% GLL points are assumed as nodes by default.
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 


pMax=2*n-1;                                                                 % Higest order of polynomials
p=0:pMax;                                                                   % orders of basis (high to low)
nP=numel(p);                                                                % number of terms in basis

etaN=inttable(n).GLL.xi;                                                    % nodal positions (GLL points)
% etaN=linspace(-1,1,nNodes);

T1=repmat(etaN',1,nP).^repmat(p,n,1);                                       % transfer matrix, functions
T2=repmat(p,n,1).*repmat(etaN',1,nP).^repmat(p-1,n,1);                      % transfer matrix, derivative
T2(:,1)=0;                                                                  % remove derivatives of constant terms
T=[T1;T2];                                                                  % transfer matrix


E=[eye(n);zeros(n)];                                                        % right hand side (conditions for shape functions)

A=T\E;                                                                      % coefficients of nodal shape functions

Psi=(eta.^p)';                                                              % values of basis functions at requested position

shp=A'*Psi;

dPsi=(p.*(eta.^(p-1)))';                                                    % values of basis functions at requested position
dPsi(1)=0;                                                                  % remove derivative of constant term
shpDeta=A'*dPsi;

shp=[shp';shpDeta'];                                                        % assemble shape functions and derivatives

