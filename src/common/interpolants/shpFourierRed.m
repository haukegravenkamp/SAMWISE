function shp = shpFourierRed(nNodes,x,~)
%% Reduced (incomplete) Fourier shape functions
% These shape functions only consist of linear terms and sine terms (no
% cosine). Hence, their convergence is suboptimal.
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
N0=1/2*(1-x);                                                               % linear functions
N1=1/2*(1+x);                                                               

n=3:nNodes;                                                                 % number of high order functions

N     = sin(pi/2*(n-2)*(x+1));                                              % Fourier shape functions 
Ndeta = cos(pi/2*(n-2)*(x+1)).* (pi/2*(n-2));                               % derivatives

shp=[N0 N N1; -0.5 Ndeta 0.5];                                              % assemble
