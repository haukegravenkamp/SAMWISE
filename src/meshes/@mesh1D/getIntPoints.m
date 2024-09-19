function [eta1D,w1D] = getIntPoints(intTab,p)

eta1D  = intTab(p+1).Gauss.xi;                                              % local coordinates at all Gauss points
w1D = intTab(p+1).Gauss.wgt;                                                % corresponding weights


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
