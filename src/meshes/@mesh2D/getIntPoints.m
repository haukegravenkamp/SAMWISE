function [eta,wgt] = getIntPoints(intTab,p)

eta1  = intTab(p(1)+1).Gauss.xi;                                            % local coordinates at all Gauss points
eta2  = intTab(p(2)+1).Gauss.xi;                                            % local coordinates at all Gauss points
wgt1  = intTab(p(1)+1).Gauss.wgt;                                           % corresponding weights
wgt2  = intTab(p(2)+1).Gauss.wgt;                                           % corresponding weights


[eta1,eta2] = meshgrid(eta1, eta2);                                         % grid of local coordinates
eta = [eta1(:), eta2(:)];
[wgt1,wgt2] = meshgrid(wgt1, wgt2);                                         % grid of weights
wgt = wgt1(:).*wgt2(:);

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
