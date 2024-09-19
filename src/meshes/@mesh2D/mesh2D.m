classdef mesh2D < mesh
    % mesh2D
    % two-dimensional mesh used for all 3D waveguides

    properties (Constant)
        nDimensions = 2                                                     % number of spatial dimensions of the mesh
    end

    properties
        edges                                                               % each row defines an edge by its vertex numbers
        edgeLengths                                                         % length of each edge
        edgeMaterials                                                       % materials of each edge
        ele2Edges                                                           % each row has the edge numbers for one element
        edges2Ele                                                           % each row has the element numbers for one edge
        edges2curve                                                         % curve number for each edge
        edgeNodes                                                           % node numbers within each edge
        nEdgeNodes                                                          % number of edge nodes for each element
        curvedEdgeFunc                                                      % function defining a curved edge, usually copied from geometry object
        curvedEdgeFuncDy                                                    % derivative of curvedEdgeFunc, computed automatically
    end

    methods

    end

    methods(Static)
        [eta1D,w1D] = getIntPoints(eta1D,w1D)
        [N,Nd1,Nd2,Na,Nad1,Nad2] = getShapeFuntions(p,eta,dof,intTab,nModes,nDom,nV)
        [N,Nd1,Nd2,Na,Nad1,Nad2] = shapeFuntionsQuad(p,eta12,dof,intTab,nModes,nDom)
        [N,Nd1,Nd2,Na,Nad1,Nad2] = shapeFuntionsTri(p,eta12,dof,intTab,nModes,nDom)
        eleMat = LoopGaussPoint(N,Ndeta1,Ndeta2,bxL,byL,bzL,brL,bxR,byR,bzR,brR,coeff,eta,xyz,nTheta,J,Jdet,wgts,nDofs,nM,matCompute)
        shpDy = getNodalDerivatives(p,dofsGl,nEleDofs,nDofs,dof,intTab,Jdet)
        [v,e] = mesh_regular2D(n1,n2,L1,L2)
        coorddomain = trianglesIntcoord(vCoord,p,v)
        femMesh = importMesh_ansys2D(fileName)
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
