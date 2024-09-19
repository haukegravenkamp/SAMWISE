function  [obj, msh] = createMesh(obj,msh,~)
%% cross-section mesh of full cylinder (y,z)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
R = obj.R;

%% vertices
theta = linspace(0,2*pi,9);                                                 % angles of vertices in 45° steps
theta = theta(1:end-1);

coord = zeros(17,3);
% inner 
coord(1:8,3)  = sin(theta)*R/2;                                             % vertical coordinate z
coord(1:8,2)  = cos(theta)*R/2;                                             % horizontal coordinate y
% outer
coord(9:16,3) = sin(theta)*R;
coord(9:16,2) = cos(theta)*R;

connectivity = [...
    17 1 2 3;...
    17 3 4 5;...
    17 5 6 7;...
    17 7 8 1;...
    1 9 10 2;...
    3 2 10 11;...
    4 3 11 12;...
    5 4 12 13;...
    6 5 13 14;...
    7 6 14 15;...
    8 7 15 16;...
    1 8 16 9
    ];

obj.curvedEdgeFunc = {@(y) sqrt(R^2-y.^2),@(y) -sqrt(R^2-y.^2)};            % define curve for all nodes as z=f(y) at R

msh.connectivity = connectivity;
msh.coord  = coord;
msh.coordV = coord;
msh.curvedEdgeFunc = obj.curvedEdgeFunc;


end
