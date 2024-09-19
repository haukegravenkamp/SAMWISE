function glb = getFEMmatrices(obj,pdes,mat,intTab,opt)

nTheta = opt.model.circumferentialOrders;

glb(numel(nTheta)) = glbSys;                                                % create one global system for each circumferential order

for i = 1: numel(nTheta)

    glb(i) = computeMatrices(obj,pdes,mat,intTab,opt,nTheta(i));            % call general matrix computation for 1D meshes

end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
