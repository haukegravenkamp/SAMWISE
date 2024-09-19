function shp = shpFourier(nNodes,x,inttable)
%% Nodal Fourier shape functions
% These shape functions consist of linear, sine, and cosine functions.
% The basis functions are converted to nodal shape functions.
% This requires the solution of a linear system of equations.
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
L=1;                                                                        % element length

s=(x+1)/2*L;                                                                % element parametrisation (0->L)

m=round((nNodes-2)/2+1e-6);                                                 % number of sine and cosine terms

nCos=m;                                                                     % number of cosine terms
nSin=m;                                                                     % number of sine terms 
if round(mod(nNodes,2))==1                                                  % if number of nodes is odd
   nCos=nCos-1;                                                             % use one sine term less
end

sn=(inttable(nNodes).GLL.xi+1)/2*L;                                         % nodal values
% sn=linspace(0,L,nNodes);

%% LINEAR FUNCTIONS
Psi1=(L-s)/L;
Psi2=s/L;

%% HARMONIC FUNCTIONS

PsiC= cos((1:nCos)'*pi*s/L);
PsiS= sin((1:nSin)'*pi*s/L);

Psi=[Psi1;PsiC;PsiS;Psi2];

T1=(L-sn')/L;                                                               % transfer matrix
Tc=cos(sn'*(1:nCos)/L*pi);
Ts=sin(sn'*(1:nSin)/L*pi);
T2=sn'/L;
T=[T1,Tc,Ts,T2];

shp=(Psi'/T);                                                               % solve system of eqs. to convert to nodal functions

%% DERIVATIVES

Psi1 = -1/L;
Psi2 =  1/L;
PsiC = -sin((1:nCos)'*pi*s/L) .* (1:nCos)'*pi/L;
PsiS =  cos((1:nSin)'*pi*s/L) .* (1:nSin)'*pi/L;
Psi=[Psi1;PsiC;PsiS;Psi2];

shpdeta=L/2*(Psi'/T);                                                       % solve system of eqs.

shp=[shp;shpdeta];                                                          % assemble shape functions and derivatives




