function [figU,ax,omegaKmode] = plotModeshape(obj,msh,omegaK,opt,mat,bcd)

% choose default options if not provided
if (nargin < 4) || ~isa(mat,'option')
    opt = option;
end
if (nargin < 5)
    mat = [];
end
if (nargin < 6) || ~isa(bcd,'bc')                                           % boundary conditions not provided (required for unbounded domains)
    bcd = bc;                                                               % set default
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% call more general plotting routine with the linePlot option
linePlot = true;
[figU,ax,omegaKmode] = plotWavefields(obj,msh,omegaK,opt,mat,bcd,linePlot);

end
