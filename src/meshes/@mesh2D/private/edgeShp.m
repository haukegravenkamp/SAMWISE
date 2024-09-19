

function [shpSide,shpSideD,modes]=edgeShp(coorSide,p,nModes,modes,intTab)

eta=coorSide*2-1;
% According to Fig. 2.24 in Provatidis, C. (2018). Precursors of Isogeometric Analysis.
% For this function, the local coordinate is scaled to -1...1

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
shpSec = shpLagrange(p+1,eta,intTab);

shpSide  = zeros(1,nModes);                                                 % shape functions
shpSideD = zeros(1,nModes);                                                 % derivatives
shpSide(modes)  = shpSec(1,:);
shpSideD(modes) = shpSec(2,:)*2;                                            % note that this is the derivative wrt the element coordinate (u,v,w)
% scaled 0...1

end
