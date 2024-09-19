function [eta,wgt]=intPoints1D(obj,nP,intType)

if nargin<2||isempty(nP)
    nP=numel(obj);                                                          % number of points, choose number of integration tables
end
if nargin<3||isempty(intType)
    intType=0;                                                              % integration type Gauss
end


if intType
    eta = obj(nP).GLL.xi;                                                   % read local coordinates
    wgt = obj(nP).GLL.wgt;                                                  % read weights
else
    eta = obj(nP).Gauss.xi;                                                 % read local coordinates
    wgt = obj(nP).Gauss.wgt;                                                % read weights
end
eta=eta(:);                                                                 % integration points columnwise
wgt=wgt(:);                                                                 % weights

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
