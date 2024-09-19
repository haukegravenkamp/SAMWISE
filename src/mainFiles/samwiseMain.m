function [geo, mat, bcd, sol, opt, res, msh, glb] = samwiseMain(geo, mat, bcd, sol, opt, res)

mat = analyzeMaterial(mat,geo);                                             % check material properties and set missing values

[geo, msh] = getMesh(geo,opt);                                              % get mesh (read from file or create)

[msh, sol.fMax] = assignElementOrder(msh, geo, mat, opt, sol.fMax);         % assign element orders

[msh,pdes] = assignPDE(msh, mat, geo, bcd);                                 % associate each element with one of the available PDEs

itg = inttable(msh.maxOrder);                                               % get integration table

msh = setDofs(msh,pdes,itg);                                                % assign degrees of freedom to each element

plotMesh(msh,opt);                                                          % plot mesh       

glb = getFEMmatrices(msh,pdes,mat,itg,opt);                                 % compute FE-matrices

glb = fluidStructureCoupling(msh,glb,pdes,mat,itg,opt);                     % add terms coupling acoustic and elastic elements

[bcd, glb] = analyzeBoundaryConditions(bcd,geo,msh,mat,opt,glb);            % set basic properties of boundary conditions

[sol,bcd,glb] = callSolver(sol,glb,msh,geo,bcd,mat,opt);                    % call solver

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
