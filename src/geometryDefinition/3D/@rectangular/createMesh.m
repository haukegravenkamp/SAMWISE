function  [obj, msh] = createMesh(obj,msh,~)

Ly = obj.Ly;
Lz = obj.Lz;


%% vertices
%    4_______ 3
%    |       |
%    |       | 
%    |       |
%    |_______|          ^ z
%    1       2          | -> y 

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
coord = [...
    0, 0,  0;...
    0, Ly, 0;...
    0, Ly, Lz;...
    0, 0,  Lz];

connectivity = [1 2 3 4];

msh.connectivity = connectivity;
msh.coord  = coord;
msh.coordV = coord;


end
