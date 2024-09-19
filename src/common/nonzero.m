function a=nonzero(a)
%% NONZERO
% Returns the entries of A that are nonzero.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
a=a(a>0);
