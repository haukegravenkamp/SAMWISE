function [eta,wgt]=intPoints2D(obj,nP1,nP2,intType)

if nargin<2||isempty(nP1)
    nP1=numel(obj);                                                         % number of points, direction 1
end
if nargin<3||isempty(nP2)
    nP2=numel(obj);                                                         % number of points, direction 2
end
if nargin<4||isempty(intType)
    intType=0;                                                              % integration type Gauss
end


if intType
    [eta1,eta2] = meshgrid(obj(nP1).GLL.xi, obj(nP2).GLL.xi);               % grid of local coordinates
    [wgt1,wgt2] = meshgrid(obj(nP1).GLL.wgt, obj(nP1).GLL.wgt);             % grid of weights
else
    [eta1,eta2] = meshgrid(obj(nP1).Gauss.xi,  obj(nP2).Gauss.xi);          % grid of local coordinates
    [wgt1,wgt2] = meshgrid(obj(nP1).Gauss.wgt, obj(nP1).Gauss.wgt);         % grid of weights
end
eta=[eta1(:),eta2(:)];                                                      % integration points columnwise

wgt=wgt1(:).*wgt2(:);                                                       % weights

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
