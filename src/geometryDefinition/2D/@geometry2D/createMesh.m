function  [obj, msh] = createMesh(obj,msh,opt)

obj = hRefine(obj,opt);                                                     % update number of elements in case of h-refinement

eleNumbers = [obj.layers.nEle];                                             % number of elements for each layer

nEle = sum(eleNumbers);                                                     % total number of elements in mesh

layerNumber = zeros(nEle,1);                                                % initialize layer numbers
materialNumber = zeros(nEle,1);                                             % initialize material numbers
orders = zeros(nEle,1);                                                     % initialize element orders
eleSize = zeros(nEle,1);                                                    % initialize element sizes

nLayers = numel(obj.layers);                                                % number of layers

eleNoBounds = cumsum(eleNumbers);                                           % accumulated sum

eleToLayer = [[1,eleNoBounds(1:end-1)+1]; eleNoBounds];                     % for each layer, first and last element number

zPosition = obj.r0;                                                         % current vertical position, updated after every layer

coord = zeros(nEle+1,3);                                                    % initialize nodal coordinates

for iLay = 1:nLayers                                                        % loop layers

    indEle = eleToLayer(1,iLay):eleToLayer(2,iLay);                         % indices of elements belonging to layer
    layerNumber(indEle) = iLay;                                             % assign layer number to elements
    materialNumber(indEle) = obj.layers(iLay).material;                     % assign material number
    eleSize(indEle) = obj.layers(iLay).thickness/obj.layers(iLay).nEle;     % assign element sizes

    layerOrder = obj.layers(iLay).eleOrder;                                 % element order requested by layer
    if isinteger(layerOrder) || isfloat(layerOrder)                         % element order is prescribed as integer or float
        orders(indEle) = round(layerOrder);                                 % use given value and round just in case
    elseif ~ischar(layerOrder) || ~strcmp(layerOrder,'auto')                % eleOrder is not 'auto' nor a float or int 
        warning(['I dont understand the eleOrder property of layer' ...
            ' no. ',num2str(iLay),'. Will try to set it automatically.'])   % display warning but continue
    end

    nodeNumbers = [indEle,(indEle(end)+1)];                                 % nodes belonging to current layer
    coord(nodeNumbers,3) = zPosition + [0;cumsum(eleSize(indEle))];

    zPosition = zPosition + obj.layers(iLay).thickness;


end

connectivity = [(1:nEle)',(2:nEle+1)'];                                     % node connectivity (including only endpoints, higher order set later)
% we assume here that all layers are connected; this is later modified 
% where neighboring elements use different PDEs

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
nNodes = nEle + 1;                                                          % preliminary number of nodes (before increasing order)

%% properties set by this function
msh.layerNumber = layerNumber;
msh.materialNumber = materialNumber;
msh.eleOrder = orders;
msh.size = eleSize;
msh.nEle = nEle;
msh.connectivity = connectivity;
msh.nNodes = nNodes;
msh.coord = coord;
msh.coordV = coord;
msh.nEleVertices = ones(nEle,1)*2;                                          % number of vertices for each element


end
