classdef mesh
    % mesh
    % general mesh class
    
    properties
        nEle                                                                % number of elements
        coord                                                               % all nodal coordinates
        coordV                                                              % vertex coordinates
        nDofs                                                               % number of degrees of freedom
        dofsLocal                                                           % local dofs of each element
        dofsGlobal                                                          % global dofs of each element
        eleOrder                                                            % order of each element
        edgeOrder                                                           % polynomial order for each edge (2D meshes)
        eleOrderAll                                                         % for 2D: p direction 1, p dir 2, p edges
        maxOrder                                                            % maximum element order
        size                                                                % size of each element (length or sqrt(surface))
        materialNumber                                                      % material assigned to each element
        pdeNumber                                                           % number of PDE for each element
        nEleNodes                                                           % number of nodes for each element
        nNodes                                                              % total number of nodes in mesh
        connectivity                                                        % node connectivity
        nEleVertices                                                        % number of vertices for each element
        pdes                                                                % storage for all PDEs
        nPde                                                                % number of PDEs
        pde2elements                                                        % for each pde stores the element numbers
        pdeDofs                                                             % list of set dofs for each PDE
        material2PDE                                                        % for each material the pde number it uses
        nDomainNodes                                                        % number of domain nodes for each element (for 2D meshes)
        elementIsCurved                                                     % indicates for each element whether it has at least one curved edge

    end
    
    methods
        function maxOrder = get.maxOrder(obj)
            maxOrder = max([obj.eleOrder(:);obj.edgeOrder(:)]);

        end

    end
    methods (Static)
        fMax = getFmax(fMax,eleSizes,cMats,materialNumber)
        [K,indb,varargout] = reduceBandwidth(K,varargin)
    end


end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
