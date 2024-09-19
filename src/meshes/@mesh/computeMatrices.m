function glb = computeMatrices(obj,pdes,mat,intTab,~,nTheta)

if nargin < 6
    nTheta = 0;                                                             % circumferential order for cylinders
end

%% interface
nPde = obj.nPde;                                                            % number of used PDEs
pdeNumber = obj.pdeNumber;                                                  % pde number of each element
glb = glbSys;                                                               % initialize global system
eleOrder = obj.eleOrder;                                                    % element order per element
eleOrderAll = obj.eleOrderAll;                                              % element order per element
materialNumber = obj.materialNumber;                                        % material number per element
connectivity = obj.connectivity;                                            % connectivity table
nDofs = obj.nDofs;                                                          % total number of dofs
dofsGlobal = obj.dofsGlobal;                                                % global dofs per element
nEleNodes = obj.nEleNodes;                                                  % number of nodes for each element
nDom = obj.nDomainNodes;
nEle = obj.nEle;
coord = obj.coord;
nEleVertices = obj.nEleVertices;
elementIsCurved = obj.elementIsCurved;
if isempty(elementIsCurved)
    elementIsCurved = false(nEle,1);
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if isempty(nDom)
    nDom = zeros(nEle,1);
end

pdeNdofs = [pdes.dof].';                                                    % number of dofs for each PDE
nCoeff = sum((nEleNodes .* pdeNdofs(pdeNumber)).^2);                        % number of nonzero entries in FEM matrices
feMatrices = glb.feMatrices;                                                % global FE matrices, presumably empty at this point
feMatrices{4,3} = [];                                                       % allocate
% assume for now only quadratic PDEs-> const, lin (2x), quad in k and 
% const, lin, quad in omega
iCoeff = zeros(nCoeff,1);                                                   % row index for sparse matrix
jCoeff = zeros(nCoeff,1);                                                   % column index for sparse matrix
mCoeff = zeros(nCoeff,numel(feMatrices));                                   % coefficients of all matrix
coeffCounter = 0;

%% loop over PDEs and compute matrices for all elements that use this PDE
% so that we don't need to redefine the element routines too often

element2Pde = (pdeNumber == (1:nPde));                                      % logical array of PDEs assigned to elements

for iPde = 1:nPde                                                           % loop PDEs

    pdeCurrent = pdes(iPde);
    dof = pdeCurrent.dof;
    eleNos = find(element2Pde(:,iPde));                                     % elements belonging to this PDE

    nTerms = pdeCurrent.nTerms;                                             % number of terms in current PDE

    bMatrices = pdeCurrent.bMatrices;                                       % read b matrices
    bxL = bMatrices(:,1);                                                   % extract matrices for x and y derivatives
    byL = bMatrices(:,2);                                                   % left (L) and right (R) refer to test and trial functions
    bzL = bMatrices(:,3);
    brL = bMatrices(:,4);
    bxR = bMatrices(:,5);
    byR = bMatrices(:,6);
    bzR = bMatrices(:,7);
    brR = bMatrices(:,8);

    ordersSpTi = pdeCurrent.ordersSpaceTime;                                % read required order of spatial and temporal derivatives

    nMatrices = (ordersSpTi(:,1)+1).*(ordersSpTi(:,2)+1);                   % number of matrices to compute for each term
    % Here, we include both first-order terms for the computation, such as
    % E12'-E11, though they will be combined into a single FE matrix.
    maxNmatrices = max(nMatrices);                                          % maxmimum number of matrices in any term
    matCompute = false(nTerms,maxNmatrices);                                % which matrices to compute for each term

    for iTerm = 1:nTerms                                                    % loop terms
        maxC = max(ordersSpTi(iTerm,1:2)+1);                                % maximum number of contributions to derivatives
        ttrial = false(maxC,1);                                             % indices of terms with derivatives in trial functions
        ttest  = false(maxC,1);                                             % indices of terms with derivatives in test  functions
        ttrial(1:ordersSpTi(iTerm,1)+1) = true;                             % we need one matrix for each spatial order less or equal requested on
        ttest(1:ordersSpTi(iTerm,2)+1) = true;                              % same for test  functions
        termCombinations = reshape(ttrial*ttest',[],1);                     % indices of combinations
        nCombinations = numel(termCombinations);                            % number of combinations
        matCompute(iTerm,1:nCombinations) = termCombinations;               % which matrices need to be computed

        % for terms without spatial derivatives, we set bx = 1, by = 0
        % so that we can use the same routine for all matrices
        if isempty(bxL{iTerm}); bxL{iTerm} = 1; end
        if isempty(byL{iTerm}); byL{iTerm} = 0; end
        if isempty(bzL{iTerm}); bzL{iTerm} = 0; end
        if isempty(brL{iTerm}); brL{iTerm} = 0; end
        if isempty(bxR{iTerm}); bxR{iTerm} = 1; end
        if isempty(byR{iTerm}); byR{iTerm} = 0; end
        if isempty(bzR{iTerm}); bzR{iTerm} = 0; end
        if isempty(brR{iTerm}); brR{iTerm} = 0; end

    end

    % LOOP ELEMENTS
    for i = 1:numel(eleNos)                                                 

        iEle = eleNos(i);                                                   % actual element number
        pAll = eleOrderAll(iEle,:);                                         % current element order
        nNodes = nEleNodes(iEle);                                           % number of nodes in element
        nV = nEleVertices(iEle);
        dofsGl = nonnan(dofsGlobal(iEle,:));
        nEleDofs = nNodes * dof;                                            % number of dofs in current element
        
        [eta,wgts] = obj.getIntPoints(intTab,eleOrder(iEle,:));             % integration points and weights depending on mesh dimension

        [N,Nd1,Nd2,Na,Nad1,Nad2] = obj.getShapeFuntions(pAll,eta,dof,intTab,nEleNodes(iEle),nDom(iEle),nV);

        % for debugging
        % plotShapefunctions(obj,pAll,nV,nDom(iEle));

        nodes = nonnan(connectivity(iEle,:));                               % node numbers
        coordNodes = coord(nodes,:);                                        % node coordinates

        [xyz,J,Jdet] = coordinateTransform(obj,nV,coordNodes,eta,N,Nd1,Nd2,~elementIsCurved(iEle));          % global coordinates and Jacobian determinant at all Gauss points
        matNo = materialNumber(iEle);                                       % material number
        coeffs = getCoefficient(pdeCurrent,mat(matNo),xyz);                 % get coefficients at all Gauss points

        [cj,ci] = meshgrid(dofsGl,dofsGl);                                  % row and column numbers in sparse matrices

        indCoeff = coeffCounter + (1:nEleDofs^2);
        iCoeff(indCoeff) = ci;
        jCoeff(indCoeff) = cj;
        coeffCounter = coeffCounter + nEleDofs^2;

        % LOOP TERMS
        for iTerm = 1:nTerms    
            
            % LOOP GAUSS POINTS
            eleMat = obj.LoopGaussPoint(Na,Nad1,Nad2,...
                bxL{iTerm},byL{iTerm},bzL{iTerm},brL{iTerm}, ...
                bxR{iTerm},byR{iTerm},bzR{iTerm},brR{iTerm},...
                coeffs{iTerm},eta,xyz,nTheta,J,Jdet,wgts,...
                nEleDofs,nMatrices(iTerm),matCompute(iTerm,:));             % get stack of matrices corresponding to current term

            for iS = 1:size(matCompute,2)
                orderTime = ordersSpTi(iTerm,3);                            % order in time
                linInd = orderTime*4+iS;
                if matCompute(iTerm,iS)
                    mCoeff(indCoeff,linInd) = reshape(eleMat(:,:,iS),[],1);
                end
            end
        end

        glb.shpDy{iEle} = obj.getNodalDerivatives(eleOrder(iEle,:),dofsGl,nEleDofs,nDofs,dof,intTab,Jdet);

    end

end

for i = 1:numel(feMatrices)
    feMatrices{i} = (sparse(iCoeff,jCoeff,mCoeff(:,i)));
end

glb.feMatrices = feMatrices;





end
