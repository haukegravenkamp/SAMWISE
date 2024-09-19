function obj=corrCurveQuad(obj,intTab)

%% CORRECTIONS FOR CURVED EDGES (QUADRILATERAL ELEMENT)
% Blending formulation as presented in
% Stipcich, G., & Piller, M. (2015). IJNME, 102, 22–43.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if isempty(obj.edges2curve)                                                 % if there are no curved edges
    return                                                                  % do nothing
end


%% interface
edges2curve = obj.edges2curve;
nEle = obj.nEle;
elementIsCurved = obj.elementIsCurved;
connectivity = obj.connectivity;
coord = obj.coord;
ele2Edges = obj.ele2Edges;
edges = obj.edges;
edgeNodes = obj.edgeNodes;
eleOrder = obj.eleOrder;
nEdgeNodes = obj.nEdgeNodes;

signF=[-1; +1; +1; -1];                                                     % signs occuring in function F

for iEle = 1:nEle
    p = eleOrder(iEle,:);
    if ~elementIsCurved(iEle) || ~any(p)>1
        continue
    end

    eleEdges = ele2Edges(iEle,:);                                           % edges belonging to current element

    nEnodesC = sum(nEdgeNodes(abs(eleEdges)));

    conC = nonnan(connectivity(iEle,:));                                  % connectivity of current element
    indVertices = conC(1:4);                                                % indices of interior nodes
    indInterior = conC((4+nEnodesC+1):end);                                 % indices of interior nodes

    vCoord = coord(indVertices,2:3);                                      % vertex coordinates

    x1 = vCoord(1,1);
    x2 = vCoord(2,1);
    x3 = vCoord(3,1);
    x4 = vCoord(4,1);
    y1 = vCoord(1,2);
    y2 = vCoord(2,2);
    y3 = vCoord(3,2);
    y4 = vCoord(4,2);

    vx = [x1 x2; x2 x3; x4 x3; x1 x4];                                      % vertex coordinates for linear interpolation
    vy = [y1 y2; y2 y3; y4 y3; y1 y4];

    eta1 = intTab(p(1)+1).GLL.xi;                                           % nodes at GLL points
    eta2 = intTab(p(2)+1).GLL.xi;

    [eta1,eta2] = meshgrid(eta1(2:end-1),eta2(2:end-1));                    % 2D grid of local coordinates
    eta1 = eta1(:);                                                         % convert to vector
    eta2 = eta2(:);


    for iEdge = 1:4                                                         % loop over edges
        edgeNr = abs(eleEdges(iEdge));                                      % current global edge number
        if ~edges2curve(edgeNr)
            continue
        end
        edgeSign = (eleEdges(iEdge)>0);                                     % sign of current edge number

        vertices = edges(edgeNr,:);                                         % vertex numbers of edge
        vCoord = coord(vertices,2:3);                                       % vertex coordinates
        edgeNodesC = nonnan(edgeNodes(edgeNr,:));                           % inner edge numbers
        eCoord = coord(edgeNodesC,2:3);                                     % inner edge coordinates
        edgeCoord=[vCoord(1,:);eCoord;vCoord(2,:)];                         % all edge node coordinates

        if ~edgeSign                                                        % edge number is negative
            edgeCoord = edgeCoord(end:-1:1,:);                              % invert order
        end
        if iEdge>2                                                          % edge 3 or 4
            edgeCoord = edgeCoord(end:-1:1,:);                              % invert order
        end
        xn = edgeCoord(:,1);                                                % nodal x coordinates
        yn = edgeCoord(:,2);                                                % nodal y coordinate

        Nx = numel(xn);                                                     % number of nodes
        for i = 1 : numel(eta1)
            if (iEdge == 1) || (iEdge == 3)
                coor1 = eta1(i);                                            % first coordinate of current point (along edge)
                coor2 = eta2(i);                                            % second coordinate of current point
            else
                coor1 = eta2(i);                                              % first coordinate of current point (along edge)
                coor2 = eta1(i);                                              % second coordinate of current point
            end

            N = shpLagrange(Nx,coor1,intTab);                               % edge shape functions
            N = N(1,:);
            s = [N*xn,N*yn];                                                % coordinates on edge

            Is=(1-coor1)/2*[vx(iEdge,1),vy(iEdge,1)]+...
                (1+coor1)/2*[vx(iEdge,2),vy(iEdge,2)];                      % linear interpolation along edge
            F=(1+signF(iEdge)*coor2)/2;                                     % linear interpolation along other direction
            c=(s-Is)*F;                                                     % correction for curvature of current edge
            coord(indInterior(i),2:3) = coord(indInterior(i),2:3)+c;              % correct node position
        end
    end
end

%% output
obj.coord = coord;

end
