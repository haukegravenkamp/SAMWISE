function [obj, pdes] = assignPDE(obj,mat,geo,bcd)
% select pde based on element and material properties

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

%% interface
pdeNumber = obj.pdeNumber;
materialNumber = obj.materialNumber;
nEle = obj.nEle;
nMat = numel(mat);                                                          % number of materials

%% check which material belongs to which pde
% get material numbers associated with boundary conditions, if any
matBc = zeros(numel(bcd),1);
for ib = 1:numel(bcd)
    if isprop(bcd(ib),'material') && ~isempty(bcd(ib).material)
        matBc(ib) = bcd(ib).material;
    end
end
materialNumber = [materialNumber;nonzero(matBc)];


pdes = pde.getAllPdes(geo);                                                 % array of all PDEs

nPdeAll = numel(pdes);                                                      % number of PDEs
material2PDE = zeros(nMat,1);                                               % which material corresponds to which PDE
pdeCounter = 0;                                                             % keep track of how many PDEs have been used
pdeUsed = zeros(nPdeAll,1);                                                 % keep track of wich PDEs have been used
pdeIndices = zeros(nPdeAll,1);                                              % indices of used PDEs
for i = 1:numel(materialNumber)                                             % loop used materials
    iMat = materialNumber(i);                                               % get element number
    for  iPDE = 1:nPdeAll                                                   % loop PDEs
        if selectPDE(pdes(iPDE),mat(iMat))                                  % check criteria of using PDE, defined by their class
            material2PDE(iMat) = iPDE;                                      % if fulfilled, assign this PDE
            if pdeUsed(iPDE) == 0                                           % if PDE hasn't been used before
                pdeCounter = pdeCounter + 1;                                % count this pde
                pdeUsed(iPDE) = pdeCounter;                                 % assign number in order of usage
                pdeIndices(pdeCounter) = iPDE;                              % store current index of used PDE
            end
        end
    end
end
pdes = pdes(nonzero(pdeIndices));                                           % store only used PDEs
nPde = numel(pdes);                                                         % update number of PDEs
material2PDECopy = material2PDE;
for iPDE = 1:nPdeAll                                                        % loop all PDEs
    material2PDE(material2PDECopy==iPDE) = pdeUsed(iPDE);                   % translate to updated number containing only used PDEs
end

if isempty(pdeNumber)                                                       % if PDE numbers per element not already defined
    pdeNumber = zeros(nEle ,1);                                             % initialize
end

for iEle = 1:nEle                                                           % loop elements

    if pdeNumber(iEle)>0                                                    % PDE number already defined (by some other means to be implemented)
        continue                                                            % do nothing
    end
    matNo = materialNumber(iEle);                                           % material number of current element
    pdeNumber(iEle) = material2PDE(matNo);                                  % read corresponding PDE number

end

element2Pde = (pdeNumber == (1:nPde));                                      % logical array of PDEs assigned to elements
pde2elements{nPde} = [];
for iPde = 1:nPde                                                           % loop PDEs
    pde2elements{iPde} = find(element2Pde(:,iPde));                         % elements belonging to this PDE
end

%% output
obj.pdes = pdes;                                                            % store PDEs with mesh object
obj.pde2elements = pde2elements;
obj.pdeNumber = pdeNumber;
obj.material2PDE = material2PDE;

end
