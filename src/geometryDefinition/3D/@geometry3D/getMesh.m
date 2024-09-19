function  [obj, msh] = getMesh(obj,opt)

msh = chooseMeshType(obj);                                                  % create mesh object

[obj, msh] = createMesh(obj,msh,opt);                                       % read or create mesh based on geometry object

msh.coord = msh.coord.*[1,obj.scaleGeometry];
msh.coordV = msh.coordV.*[1,obj.scaleGeometry];


nEle = size(msh.connectivity,1);                                            % total number of elements in mesh
nNodes = size(msh.coord,1);
msh.nNodes = nNodes;
msh.nEle = nEle;

msh.nEleVertices = sum(msh.connectivity>0,2);                               % number of vertices for each element

%% element materials
materialNumber = obj.materialNumber;                                        % material numbers as provided by geometry object

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if numel(materialNumber)~=nEle                                              % not one material per element given
    if numel(materialNumber) > 1                                            % several materials given, inconsistent
        warning(['The material number(s) of the geometry should be given ' ...
            'either as single integer or a vector of material numbers of each' ...
            'element. Using the first material for the entire domain.'])
    end
    materialNumber = ones(nEle,1)*materialNumber(1);                        % use only first material
end
msh.materialNumber = materialNumber;


%% element sizes
msh = getEdgeProperties(msh);
msh.size = max(msh.edgeLengths,2);

%% refine mesh
% to do
obj = hRefine(obj,opt);                                                     % update number of elements in case of h-refinement


end
