function  [obj, msh] = createMesh(obj,msh,opt)

Ly = obj.Ly;
Lz = obj.Lz;
hy = obj.hy;
hz = obj.hz;


%% vertices
%    7____6_5
%    |____  |
%    8    |4|3
%         | |
%         |_|              ^ z
%         1 2              | -> y 
% 
          
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
coord = [...
    0, Ly-hy, 0    ;...
    0, Ly,    0    ;...
    0, Ly,    Lz-hz;...
    0, Ly-hy, Lz-hz;...
    0, Ly,    Lz   ;...
    0, Ly-hy, Lz   ;...
    0, 0,     Lz   ;...
    0, 0,     Lz-hz
    ];


connectivity = [1 2 3 4; 4 3 5 6; 8 4 6 7];

msh.connectivity = connectivity;
msh.coord  = coord;
msh.coordV = coord;


end
