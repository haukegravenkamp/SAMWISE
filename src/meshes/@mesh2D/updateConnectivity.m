function obj = updateConnectivity(obj,intTab)

%% interface
nEle             = obj.nEle;                                                % read number of elements
connectivity     = obj.connectivity;
nNodes           = obj.nNodes;
nEleNodes        = obj.nEleNodes;
eleOrder         = obj.eleOrder;
edgeOrder        = obj.edgeOrder;
ele2Edges        = obj.ele2Edges;
nEdges           = size(obj.edges,1);
nDomainNodes     = obj.nDomainNodes;
coordV           = obj.coordV;
edges            = obj.edges;
edges2curve      = obj.edges2curve;
curvedEdgeFunc   = obj.curvedEdgeFunc;
curvedEdgeFuncDy = obj.curvedEdgeFuncDy;
nEleVertices     = obj.nEleVertices;

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

%% update connectivity with interior nodes
maxNodesPerElement = max(nEleNodes);                                        % maximum number of nodes per element

coord = zeros(maxNodesPerElement*nEle,3);                                   % extend coordinate table
coord(1:nNodes,:)=coordV;

if size(ele2Edges,2) == 3                                                   % if all triangles
    maxConnect = maxNodesPerElement + 1; 
else
    maxConnect = maxNodesPerElement;
end
if size(connectivity,2)<maxConnect                                          % extend connectivity table to account for edge and domain nodes
    connectivity(end,maxConnect) = 0;
    connectivity(connectivity==0) = nan;
end

nodeCounter = nNodes;                                                       % keep track of node number

% interior nodes on edges
edgeNodes = nan(nEdges,max(edgeOrder)-1);
nEdgeNodes = zeros(nEdges,1);                                               % number of inner nodes on each edge
for iEdge = 1:nEdges                                                        % loop edges
    nEnodes = edgeOrder(iEdge)-1;                                           % number of edge nodes, not including vertices
    nEdgeNodes(iEdge) = nEnodes;
    if nEnodes == 0
        continue
    end
    
    edgeNodes(iEdge,1:nEnodes) = nodeCounter + (1:nEnodes);                 % assign new nodes to edges
    vEdge = coordV(edges(iEdge,:),2:3);                                     % edge vertices (y,z)
    y1 = vEdge(1,1);
    y2 = vEdge(2,1);
    Ly = vEdge(2,1) - vEdge(1,1);                                           % linear edge length along y
    Lz = vEdge(2,2) - vEdge(1,2);                                           % linear edge length along z
    eta = intTab(edgeOrder(iEdge)+1).GLL.xi;                                % local coordinates
    eta = eta(2:end-1);                                                     % only inner nodes
    if ~isempty(edges2curve) && edges2curve(iEdge)>0                        % if edge is curved
        fzofy = curvedEdgeFunc{edges2curve(iEdge)};                         % read function describing curvature
        df = curvedEdgeFuncDy{edges2curve(iEdge)};                          % read derivative
        [xi,yi] = curveNodePos(fzofy,df,y1,y2,eta);                         % find global coordinates on curve
        coordEdge = [xi; yi].';                                             % concatenate
    else                                                                    % no curvature
        midC = sum(vEdge)/2;                                                % edge center
        coordEdge = [midC(1)+eta*Ly/2;midC(2)+eta*Lz/2].';                  % nodes along straight edge
    end
    coord(nodeCounter + (1:nEnodes), 2:3) = coordEdge;                      % store in coordinate table
    nodeCounter = nodeCounter + nEnodes;                                    % update counter
end

% interior nodes 
for iEle = 1:nEle                                                           % loop elements
    nInteriorNodes = nDomainNodes(iEle);                                    % number of interior nodes in current element

    % if nInteriorNodes == 0                                                  % if no interior nodes present
    %     continue                                                            % do nothing
    % end
    edgesC = nonnan(ele2Edges(iEle,:));                                     % edges of current element
    nEdgesC = numel(edgesC);                                                % number of edges
    counterLocal = 4;                                                       % counter for nodes in current element, always start at 4

    for iEdge = 1:nEdgesC                                                   % loop edges of current element
        edgeNodesC = nonnan(edgeNodes(abs(edgesC(iEdge)),:));               % node numbers of current edge
        nEdgeNodesC = numel(edgeNodesC);                                    % number of nodes on edge
        ind = counterLocal + (1:nEdgeNodesC);                               % indices in connectivity of current element
        if edgesC(iEdge) < 0                                                % negative orientation of edge
            edgeNodesC = edgeNodesC(end:-1:1);                              % invert order
        end
        connectivity(iEle,ind) = edgeNodesC;                                % assign nodes for edges
        counterLocal = counterLocal + nEdgeNodesC;                          % update local counter
    end

    ind = counterLocal + (1:nInteriorNodes);                                % indices in connectivity table
    connectivity(iEle,ind) = (1:nInteriorNodes) + nodeCounter;              % assign nodes for interior

    e1 = intTab(eleOrder(iEle,1)+1).GLL.xi;                                 % local coordinate, direction 1
    coordEle = coordV(nonnan(connectivity(iEle,1:4)),:);                    % vertex coordinates

    if nEleVertices(iEle)==3                                                      % triangles
        coordGl = obj.trianglesIntcoord(coordEle(:,2:3),eleOrder(iEle,1),e1);
        coord(nodeCounter + (1:nInteriorNodes), 2:3) = coordGl;
    else                                                                    % quads
        e2 = intTab(eleOrder(iEle,2)+1).GLL.xi;                             % local coordinate, direction 2
        [eta1,eta2] = meshgrid(e1(2:end-1),e2(2:end-1));                    % interior nodes only
        eta12 = [eta1(:),eta2(:)];
        coordGl = coordinateTransform(obj,4,coordEle,eta12);                  % linear coordinate transform
        coord(nodeCounter + (1:nInteriorNodes), :) = coordGl;
    end


    nodeCounter = nodeCounter + nInteriorNodes;
end
nNodes = nodeCounter;                                                       % update total number of nodes including interior nodes

coord((nodeCounter+1):end,:) = [];

%% output
obj.connectivity = connectivity;
obj.nNodes = nNodes;
obj.coord = coord;
obj.edgeNodes = edgeNodes;
obj.nEdgeNodes = nEdgeNodes;

%% correction for elements with curved edges
obj=corrCurveQuad(obj,intTab);

end
