classdef mesh1D < mesh
    % mesh1D
    % one-dimensional mesh used for plates and cylinders

    properties (Constant)
        nDimensions = 1                                                     % number of spatial dimensions of the mesh
    end

    properties
        layerNumber                                                         % which layer the element belongs to
    end

    methods(Static)
        [eta1D,w1D] = getIntPoints(intTab,p)
        [N,Nd1,Nd2,Na,Nad1,Nad2] = getShapeFuntions(p,eta,dof,intTab,~,~,~)
        shpDy = getNodalDerivatives(p,dofsGl,nEleDofs,nDofs,dof,intTab,Jdet)
        eleMat = LoopGaussPoint(N,Ndeta1,Ndeta2,bxL,byL,bzL,brL,bxR,byR,bzR,brR,coeff,eta,xyz,nTheta,J,Jdet,wgts,nDofs,nM,matCompute)

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
