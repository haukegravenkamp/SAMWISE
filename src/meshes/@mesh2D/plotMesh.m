function h = plotMesh(obj,opt)

if nargin < 2
    opt = option;
end

if ~opt.plotting.plotMesh
    return
end

plotNodes = ~opt.plotting.plotNodes;

%% interface
coord = obj.coord;
coordV = obj.coordV;
edges = obj.edges;
curvedEdgeFunc = obj.curvedEdgeFunc;
edges2curve = obj.edges2curve;
maxOrder = obj.maxOrder;

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% plot straight edges
nEdges=size(edges,1);
ycoord=nan(2,nEdges);
zcoord=nan(2,nEdges);
for i = 1:nEdges
    if ~isempty(edges2curve) && edges2curve(i)>0
        continue
    end
   ycoord(:,i) = coordV(edges(i,:),2);
   zcoord(:,i) = coordV(edges(i,:),3);
end
h = plot(ycoord,zcoord,'-k','Linewidth',1);
hold all

%% plot curved edges
N = 50*maxOrder;
nEcurved = sum(edges2curve>0);
ycoordC=nan(N,nEcurved);
zcoordC=nan(N,nEcurved);
curveNumbers = unique(edges2curve(edges2curve>0));
edgeCount = 0;
for i = 1:numel(curveNumbers)
    iC = curveNumbers(i);
    fy = curvedEdgeFunc{iC};
    edgeInd = find(edges2curve==iC);
    for j = 1:numel(edgeInd)
        iEdge =edgeInd(j);
        y1 = coordV(edges(iEdge,1),2);
        y2 = coordV(edges(iEdge,2),2);
        y = linspace(y1,y2,N);
        z = fy(y);
        edgeCount = edgeCount +1;
        ycoordC(:,edgeCount) = y;
        zcoordC(:,edgeCount) = z;

    end

end

h = plot(ycoordC,zcoordC,'-k','Linewidth',1);
hold all
mystdfig(12,[],2)
axis equal
xlabel('$y$')
ylabel('$z$')

%% plot nodes
if ~isscalar(plotNodes)
    if size(coord,1)<500
        plotNodes = true;
    else
        plotNodes = false;
    end
end

if plotNodes
    plot(coord(:,2),coord(:,3),'.k','MarkerSize',14)
end

drawnow


end
