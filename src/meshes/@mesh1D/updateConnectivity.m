function obj = updateConnectivity(obj,intTab)

%% interface
nEle = obj.nEle;                                                            % read number of elements
connectivity = obj.connectivity;
nNodes = obj.nNodes;
nEleNodes = obj.nEleNodes;
eleOrder = obj.eleOrder;
coordV = obj.coordV;

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% update connectivity with interior nodes
maxNodesPerElement = max(nEleNodes);                                        % maximum number of nodes per element
coord = zeros(maxNodesPerElement*nEle,3);                                   % extend coordinate table
coord(1:nNodes,:)=coordV;

maxNodesPerElement = max(nEleNodes);                                        % maximum number of nodes per element
if size(connectivity,2)<maxNodesPerElement                                  % extend connectivity table to account for interior nodes
    connectivity(end,maxNodesPerElement) = 0;
    connectivity(connectivity==0) = nan;
end
nodeCounter = nNodes;                                                       % keep track of node number
for iEle = 1:nEle                                                           % loop elements
    nInteriorNodes = nEleNodes(iEle)-2;                                     % number of interior nodes in current element
    connectivity(iEle,3:nInteriorNodes+2) = ...
        (1:nInteriorNodes) + nodeCounter;                                   % assign nodes for interior

    % get coordinates
    vEdge = coordV(connectivity(iEle,1:2),2:3);                             % vertices (y,z)
    Ly = vEdge(2,1) - vEdge(1,1);                                           % linear edge length along y
    Lz = vEdge(2,2) - vEdge(1,2);                                           % linear edge length along z
    eta = intTab(eleOrder(iEle)+1).GLL.xi;                                  % local coordinates
    eta = eta(2:end-1);                                                     % only inner nodes
    midC = sum(vEdge)/2;                                                    % edge center
    coordEdge = [midC(1)+eta*Ly/2;midC(2)+eta*Lz/2].';                      % nodes along straight edge
    coord(nodeCounter + (1:nInteriorNodes), 2:3) = coordEdge;               % store in coordinate table
    nodeCounter = nodeCounter + nInteriorNodes;

end
coord = coord(1:nodeCounter,:);
nNodes = sum(eleOrder)+1;                                                   % update total number of nodes including interior nodes

%% output
obj.connectivity = connectivity;
obj.nNodes = nNodes;
obj.coord = coord;
end
