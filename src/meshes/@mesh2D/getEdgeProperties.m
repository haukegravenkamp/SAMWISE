function obj = getEdgeProperties(obj)

%% interface
coord          = obj.coord;                                                 % nodal coordinates
connect        = obj.connectivity;                                          % connectivity table
materialNumber = obj.materialNumber;                                        % material number of each element
curvedEdgeFunc = obj.curvedEdgeFunc;                                        % function handles defining curves
nEle           = obj.nEle;                                                  % number of elements

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% get edges
% assumes we have no more than four edges per element
connectCycl = [connect,connect(:,1)];                                       % vertex connectivity, cylcic
nSidesMax = size(connect,2);                                                % maximum number of sides of any element
nEdgesMax = nEle*nSidesMax;                                                 % maximum number of edges by counting shared edges twice
edges = zeros(nEdgesMax,2);                                                 % edges defined by vertices
countE = 0;                                                                 % counter
for i = 1:nEle                                                              % loop elements       
    con = nonnan(connectCycl(i,:));                                      
    for j = 1:numel(con)-1                                                  % loop edges
        countE = countE + 1;                                                % increase conuter
        edges(countE,:) = sort(con(j:j+1));                                 % edges, increasing vertex numbers by definition
    end
end
edges = edges(any(edges,2),:);                                              % throw out empty edges
edges = unique(edges,"rows");                                               % throw out duplicate edges
nEdges = size(edges,1);                                                     % number of edges

%% assign edges to elements
ele2Edges = nan(nEle,nSidesMax);                                                    % each row contains the edge numbers of an element
for i = 1:nEle                                                              % loop elements
    con = nonnan(connectCycl(i,:));                                         % connectivity of current element
    for j = 1:numel(con)-1                                                  % loop edges of current element
        [conS,indSort] = sort(con(j:j+1));                                  % connectivity of current element, increasing order
        edgeNo = find((edges(:,1) == conS(1)) & (edges(:,2) == conS(2)));   % edge containing the same vertices
        if indSort(1) == 1
            signEdge = 1;
        else
            signEdge = -1;
        end
        ele2Edges(i,j) = edgeNo*signEdge;
    end
end

%% assign elements to edges
edges2Ele = nan(nEdges,2);
for i = 1:nEdges 
    ele = find((sum((connect == edges(i,1)) + (connect == edges(i,2)),2) == 2));
    edges2Ele(i,1:numel(ele)) = ele;
end

%% assing curves
nCurves = numel(curvedEdgeFunc);                                            % number of explicitly defined curves
if nCurves > 0                                                              % if there are curves
    syms y                                                                  % symbolic variable for finding the derivative
    yC = coord(:,2);                                                        % read y coordinates
    zC = coord(:,3);                                                        % read z coordinates
    onCurveTh = 1e-6/max(abs(coord(:)));                                    % threshold for point being on the curve
    if ~iscell(curvedEdgeFunc)                                              % curve function not inside a cell (may happen if there's only one)
        curvedEdgeFunc = {curvedEdgeFunc};                                  % create cell for consistency
    end
    curvedEdgeFuncDy{nCurves} = [];                                         % allocate derivatives (needed for arc length)
    edges2curve = zeros(size(edges,1),1);                                   % allocate curve number for each edge
    for iC = 1:nCurves                                                      % loop curves
        fy = curvedEdgeFunc{iC};                                            % read function handle                              
        indCoord = abs(zC-fy(yC))<onCurveTh;                                % indices of vertices on curve
        edgesOnCurve = sum(ismember(edges,find(indCoord)),2)==2;            % edges whose vertices are both on curve
        edges2curve(edgesOnCurve) = iC;                                     % assign curve number to edge
        curvedEdgeFuncDy{iC} = matlabFunction(diff(fy(y)));                 % compute derivative
    end
else                                                                        % no curves defined
    edges2curve = [];
    curvedEdgeFuncDy = [];
end
elementIsCurved = false(nEle,1);                                            % flags indicating for each edge whether it's curved
elementIsCurved(unique(nonnan(edges2Ele(edges2curve>0,:)))) = true;


%% compute edge lengths
edgeLengths = zeros(nEdges,1);
for i = 1:nEdges
    edgeCoord = coord(edges(i,:),2:3);
    if ~isempty(edges2curve) && edges2curve(i)>0
        y1 = edgeCoord(1,1);
        y2 = edgeCoord(2,1);
        df = curvedEdgeFuncDy{edges2curve(i)};
        integrand=@(y) sqrt(1 + abs(df(y)).^2);                                     % integrand (function handle)
        edgeLengths(i) = abs(integral( integrand,y1,y2));                                           % arc length (integral from x1 to variable b)
    else
        edgeLengths(i)=sqrt(sum((edgeCoord(2,:)-edgeCoord(1,:)).^2));
    end
end

%% get materials of each edge
if isscalar(unique(materialNumber))
    edgeMaterials = ones(nEdges,1)*materialNumber(1);
else
    edgeMaterials = nan(nEdges,2);
    for i = 1:nEdges
        eleMat = unique(materialNumber(nonnan(edges2Ele(i,:))));
        edgeMaterials(i,1:numel(eleMat)) = eleMat;
    end
end

%% output
obj.edges            = edges;
obj.edgeLengths      = edgeLengths;
obj.edges2Ele        = edges2Ele;
obj.ele2Edges        = ele2Edges;
obj.edgeMaterials    = edgeMaterials;
obj.edges2curve      = edges2curve;
obj.elementIsCurved  = elementIsCurved;
obj.curvedEdgeFuncDy = curvedEdgeFuncDy;

end
