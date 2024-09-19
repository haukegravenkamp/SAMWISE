function [obj, fMax] = assignElementOrder(obj,~,mat,opt,fMax)
% assign interpolation order to elements
% affects only elements that do not have a valid element order already
% defined

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
eleOrderDefault  = opt.numerics.eleOrderDefault;
eleOrderFunction = opt.numerics.eleOrderFunction;
nEle             = obj.nEle;                                                 % number of elements
materialNumber   = obj.materialNumber;
edgeLengths      = obj.edgeLengths;
edgeMaterials    = obj.edgeMaterials;
eleSizes         = obj.size;
usedMaterials    = obj.materialNumber;
ele2Edges        = obj.ele2Edges;

%% get element orders 

% get wave velocity relevant for computing dimensionless frequency
cMats = getWaveVelocity(mat,usedMaterials);
fMax = obj.getFmax(fMax,eleSizes,cMats,materialNumber);
% get required order for each edge
edgeOrder = getEleOrder(obj,[],cMats,edgeMaterials,edgeLengths,fMax,eleOrderDefault,eleOrderFunction);

eleEdgeLength = edgeLengths(abs(ele2Edges));                                % lenghts of all edges of all elements
if size(eleEdgeLength,2)==1
    eleEdgeLength = eleEdgeLength.';
end
maxLengthdir1 = max(eleEdgeLength(:,[1,3]),[],2);                           % maximum edge length in 1st direction
if size(eleEdgeLength,2)>3                                                  % if it's not all triangles
    maxLengthdir2 = max(eleEdgeLength(:,[2,4]),[],2);                       % maximum edge length in 2nd direction
else                                                                        % all triangles
    maxLengthdir2 = maxLengthdir1;                                          % just copy lenghts of first direction
end
% get required element order for each direction
eleOrderDir1 = getEleOrder(obj,[],cMats,materialNumber,maxLengthdir1,fMax,eleOrderDefault,eleOrderFunction);
eleOrderDir2 = getEleOrder(obj,[],cMats,materialNumber,maxLengthdir2,fMax,eleOrderDefault,eleOrderFunction);

orderKnown = obj.eleOrder;                                                  % element orders already defined for some reason
nKnown = size(orderKnown,1);                                              
if nKnown<nEle                                                              % not defined for all elements
    eleOrder = zeros(nEle,2);                                               % allocate element orders
    if nKnown > 0                                                           % if order known for some elements
        eleOrder(nKnown,1:size(orderKnown,2)) = orderKnown;                 % use known values
    end
    ind = (nKnown+1):nEle;                                                  % indices of unknown elements
    eleOrder(ind,:) = [eleOrderDir1(ind),eleOrderDir2(ind)];                % assign computed values
end

nEleNodes = zeros(nEle,1);
nDomainNodes = zeros(nEle,1);
for iEle = 1:nEle
    edgeNos = unique(abs(nonnan(ele2Edges(iEle,:))));
    nEdgeNodes = sum(edgeOrder(edgeNos));
    n = eleOrder(iEle,:)-1;                                                 % number of nodes per direction
    if numel(edgeNos) == 3                                                  % triangle
        nDomainNodes(iEle) = (n(1)-1)/2*n(1);
    elseif numel(edgeNos) == 4                                              % rectangle
        nDomainNodes(iEle) = n(1)*n(2);
    end
    nEleNodes(iEle) = nEdgeNodes + nDomainNodes(iEle);
end

if nEle == 1
    edgeOrder = edgeOrder.';
end
%% output
obj.eleOrder     = eleOrder;
obj.eleOrderAll  = [eleOrder,edgeOrder(abs(ele2Edges))];
obj.nEleNodes    = nEleNodes;
obj.edgeOrder    = edgeOrder;
obj.nDomainNodes = nDomainNodes;


end
