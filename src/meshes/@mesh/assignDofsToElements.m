
function obj = assignDofsToElements(obj,pdes)


%% interface
pdeDofsPerNode = [pdes.dof];                                                % number of dofs per node in each PDE

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
eleNodalDofs = pdeDofsPerNode(obj.pdeNumber);                               % number of dofs per node of each element
eleNodalDofs = eleNodalDofs(:);                                             % make sure it's a column vector

connectivity = obj.connectivity;
pdeNumber = obj.pdeNumber;
nEleNodes = obj.nEleNodes;

nEleDofs = obj.nEleNodes(:).*eleNodalDofs;                                  % total number of dofs per element

nEle = obj.nEle;                                                            % read number of elements
nNodes = obj.nNodes; 

%% create arrays of dofs for each PDE
nPde = numel(pdes);                                                         % number of PDEs
pdeDofs = {};
for iPde = 1:nPde
    pdeDofs{iPde} = nan(nNodes,pdeDofsPerNode(iPde));
end


%% assign dofs
maxEleDofs = max(nEleDofs);                                                 % maximum number of dofs per element

dofsLocal = nan(nEle,maxEleDofs);                                       % allocate array of local and global dofs
dofsGlobal = nan(nEle,maxEleDofs);

dofsAssigned = 0;
for iEle = 1:nEle
    iPde = pdeNumber(iEle);                                                 % read pde number
    nDofNode = pdeDofsPerNode(iPde);                                        % number of dofs per node
    eleNodes = nonzero(connectivity(iEle,:));                               % read nodes belonging to current element
    nDof = nDofNode * nEleNodes(iEle);                                      % total number of dofs in element
    eleGlbDofs = zeros(nDofNode,nEleNodes(iEle));                           % allocate global dofs

    for iNode=1:numel(eleNodes)                                             % loop nodes
        dofsInPde = nonnan(pdeDofs{iPde}(eleNodes(iNode),:));               % read dofs already assigned to pde at the node
        if isempty(dofsInPde)                                               % if non assigned
            dofsToSet = (1:nDofNode) + dofsAssigned;                        % create new dofs
            pdeDofs{iPde}(eleNodes(iNode),:) = dofsToSet;                   % store with pde
            dofsAssigned = dofsAssigned + nDofNode;                         % update dof counter
        else                                                                % dofs exist already
            dofsToSet = dofsInPde;                                          % use existing dofs
        end
        eleGlbDofs(:,iNode) = dofsToSet;                                    % store dofs
    end
    % eleGlbDofs =  eleGlbDofs(:,[1,3:numel(eleNodes),2]);
    dofsGlobal(iEle,1:nDof) = eleGlbDofs(:);
    dofsLocal(iEle,1:nDof) = 1:nDof;

end

%% output

obj.dofsGlobal = dofsGlobal;
obj.dofsLocal  = dofsLocal;
obj.nDofs      = dofsAssigned;
obj.pdeDofs    = pdeDofs;
obj.nPde       = nPde;

end
