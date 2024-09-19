function Z=Zmatrix(omega,E0,E11,E12,E2,M0,C0)
%% Zmatrix
% Z-matrix, used in the notation of the SBFEM to linearize the second-order
% differential equation

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if (nargin < 6) || isempty(M0)
    M0 = 0;
end

if (nargin < 7) || isempty(C0)
    C0 = 0;
end

E0inv = E0\eye(size(E0));
Z1 =  E0inv*E12;
Z2 = -E0inv;
Z3 =  1i*omega*C0 + omega^2*M0-E2 + E11*E0inv*E12;
Z4 = -E11*E0inv;

Z = [Z1,Z2;Z3,Z4];
