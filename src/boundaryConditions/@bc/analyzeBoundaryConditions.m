function [obj,glb] = analyzeBoundaryConditions(obj,geo,msh,mat,opt,glb)
% check assigned boundary conditions are obtain relevant dofs, materials,
% etc.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
for iBc = 1:numel(obj)                                                      % loop boundary conditions
    obj(iBc) = getBcLocation(msh,obj(iBc));                                 % read mesh information depending on mesh type
end

for i = 1:numel(glb)
    glb(i) = findFixedDofs(obj,msh,glb(i));
end
end
