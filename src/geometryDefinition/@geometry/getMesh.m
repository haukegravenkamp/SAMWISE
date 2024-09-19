function  [obj, msh] = getMesh(obj,opt)

msh = chooseMeshType(obj);

[obj, msh] = createMesh(obj,msh,opt);

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
