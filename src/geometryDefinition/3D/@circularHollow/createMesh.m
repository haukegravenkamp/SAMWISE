function  [obj, msh] = createMesh(obj,msh,~)
%% cross-section mesh of full cylinder (y,z)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
Ri = obj.Ri;  % inner radius
Ro = obj.Ro;  % outer radius

%% vertices

% determine number of layers such that the element radius does not increase
% by more than a factor of 1.5 in each element
nLayer = ceil(log(Ro/Ri)/log(3/2));

[v,conn]=mesh2D.mesh_regular2D(8,nLayer,2*pi,Ro-Ri);                        % create regular mesh in (theta,r)
th = v(:,1);                                                                % angle
r = v(:,2) + Ri;                                                            % radius

indFirst = 1:nLayer+1;                                                      % indices of nodes at theta = 0
indDupl = numel(r)-nLayer:numel(r);                                         % indices of nodes at theta = 2pi
for i = 1:numel(indDupl)
    conn(conn==indDupl(i)) = indFirst(i);
end
r(indDupl,:) = [];
th(indDupl,:) = [];

coord = zeros(numel(r),3);
coord(:,2)  = cos(th).*r;                                                   % horizontal coordinate y
coord(:,3)  = sin(th).*r;                                                   % vertical coordinate z

ru = unique(r);                                                             % unique values of radius

curvedEdgeFunc{numel(ru)*2} = [];

for i = 1:numel(ru)
    curvedEdgeFunc{2*(i-1) + 1} = @(y)  sqrt(ru(i)^2-y.^2);
    curvedEdgeFunc{2*(i-1) + 2} = @(y) -sqrt(ru(i)^2-y.^2);
end

obj.curvedEdgeFunc = curvedEdgeFunc;

msh.connectivity   = conn;
msh.coord          = coord;
msh.coordV         = coord;
msh.curvedEdgeFunc = obj.curvedEdgeFunc;


end
