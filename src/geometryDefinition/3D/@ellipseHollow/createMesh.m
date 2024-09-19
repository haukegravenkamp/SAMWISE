function  [obj, msh] = createMesh(obj,msh,~)
%% cross-section mesh of full cylinder (y,z)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
Ai = obj.Ai;  % inner half width
Ao = obj.Ao;  % outer half width
Bi = obj.Bi;  % inner half height
Bo = obj.Bo;  % outer half height

%% vertices

% determine number of layers such that the element radius does not increase
% by more than a factor of 1.5 in each element (adapted from cylinder, not precise for ellipse)
nLayer = ceil(log(Ao/Ai)/log(3/2));

[v,conn]=mesh2D.mesh_regular2D(8,nLayer,2*pi,Ao-Ai);                        % create regular mesh in (theta,a)
th = v(:,1);                                                                % angle
a  = v(:,2) + Ai;                                                           % width
b  = v(:,2)/(Ao-Ai)*(Bo-Bi) + Bi;                                           % height


indFirst = 1:nLayer+1;                                                      % indices of nodes at theta = 0
indDupl = numel(a)-nLayer:numel(a);                                         % indices of nodes at theta = 2pi
for i = 1:numel(indDupl)
    conn(conn==indDupl(i)) = indFirst(i);
end
a(indDupl,:) = [];
b(indDupl,:) = [];
th(indDupl,:) = [];


coord = zeros(numel(a),3);
coord(:,2)  = cos(th).*a;                                                   % horizontal coordinate y
coord(:,3)  = sin(th).*b;                                                   % vertical coordinate z

au = unique(a);                                                             % unique values of radius
bu = unique(b);                                                             % unique values of radius


curvedEdgeFunc{numel(au)*2} = [];

for i = 1:numel(au)
    curvedEdgeFunc{2*(i-1) + 1} = @(y)  sqrt(1-y.^2/au(i)^2)*bu(i);
    curvedEdgeFunc{2*(i-1) + 2} = @(y) -sqrt(1-y.^2/au(i)^2)*bu(i);
end

obj.curvedEdgeFunc = curvedEdgeFunc;

msh.connectivity   = conn;
msh.coord          = coord;
msh.coordV         = coord;
msh.curvedEdgeFunc = obj.curvedEdgeFunc;


end
