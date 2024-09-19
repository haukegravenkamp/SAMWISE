function shp = shpFourierPure(nNodes,x,~)
%% "Pure" Fourier shape functions. 
% These shape functions consist of sine and cosine functions.
% All cosine terms are nonzero at the element end points.
% Hence, elements of this type can not be coupled straightforwardly.
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
L=2;                                                                        % element length

s=(x+1)/2*L;                                                                % element parametrisation (0->L)

m=round((nNodes-2)/2+1e-6);                                                 % number of sine and cosine terms

nCos=m;                                                                     % number of cosine terms
nSin=m;                                                                     % number of sine terms
if round(mod(nNodes,2))==1                                                  % if number of nodes is odd
    nSin=nSin-1;                                                             % use one sine term less
end

%% LINEAR FUNCTIONS

Psi1=(L-s)/L;
Psi2=s/L;

%% HARMONIC FUNCTIONS

PsiC= cos((1:nCos)'*pi*s/L);
PsiS= sin((1:nSin)'*pi*s/L);

Psi=[Psi1;PsiC;PsiS;Psi2];

shp=Psi';

%% DERIVATIVES

Psi1 = -1/L;
Psi2 =  1/L;
PsiC = -sin((1:nCos)'*pi*s/L) .* (1:nCos)'*pi/L;
PsiS =  cos((1:nSin)'*pi*s/L) .* (1:nSin)'*pi/L;
Psi=[Psi1;PsiC;PsiS;Psi2];

shpdeta=L/2*(Psi');
shp=[shp;shpdeta];




