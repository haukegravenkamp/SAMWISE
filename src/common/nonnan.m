function a=nonnan(a)
%% NONNAN
% Returns the entries of A that are not NANs.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
a=a(~isnan(a));
