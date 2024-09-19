function  [obj, msh] = createMesh(obj,msh,~)
%% cross-section mesh of full cylinder (y,z)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
fileName = obj.fileName;  % name of mesh file

femMesh = mesh2D.importMesh_ansys2D(fileName);
coordyz = femMesh.vertices;
coord = [zeros(size(coordyz,1),1),coordyz];

msh.connectivity   = femMesh.connectivity;
msh.coord          = coord;
msh.coordV         = coord;

end
