function [coordV,conne,dofsGl,matNo,sizes,isUnbounded,pde2ele,unbBTa]=...
    includeUnboundedDomain(bcd,nEle,coordV,conne,dofsGl,matNo,mat2PDE,nPde,pde2ele,sizes,leUnb,unboundedModel)

if strcmp(unboundedModel,'exact')                                           % exact incorporation of embedded boundary conditions
    nNodes = size(coordV,1);                                                % number of nodes in waveguide
    [zmin, nodeMin] = min(coordV(:,3));                                     % min z-coordinate of waveguide
    [zmax, nodeMax] = max(coordV(:,3));                                     % max z-coordinate
    coordVunb = nan(2,3);                                                   % allocate location of unbounded domain
    dofsUn = nan(2,size(dofsGl,2));                                         % allocate dofs of unbounded domains
    connUnb = nan(2,size(conne,2));                                         % allocate connectivity of unbounded domains
    matNoUnb = nan(2,1);                                                    % allocate material numbers
    unbPde = nan(2,1);                                                      % allocate PDE number of unbounded domains
    unbBT  = nan(2,1);                                                      % allocate flag for bottom and top surface
    iNew = 0;                                                               % keep track of additional elements for unbounded
    for i = 1:numel(bcd)                                                    % loop boundary conditions
        if isa(bcd(i),'bcUnbounded')                                        % unbounded
            iNew = iNew + 1;                                                % add element
            unbPde(iNew) = mat2PDE(bcd(i).material);                        % store PDE
            matNoUnb(iNew) = bcd(i).material;                               % store material number
            dofsUnbC = bcd(i).dofsUnbounded;                                % dofs of current unbounded domain
            dofsUn(iNew,1:numel(dofsUnbC)) = dofsUnbC;                      % store coupld dofs in waveguide
            if strcmp(bcd(i).location,'bottom')                             % lower surface
                connUnb(iNew,1:2) = [iNew+nNodes, nodeMin];                 % connectivity: new node below, lower surface
                unbBT(iNew) = 1;                                            % flag
                coordVunb(iNew,:) = [0 0 zmin - leUnb];                     
            elseif strcmp(bcd(i).location,'top')                            % upper surface
                connUnb(iNew,1:2) = [nodeMax, iNew+nNodes];                 % connectivity: upper surface, new node above
                unbBT(iNew) = 2;
                coordVunb(iNew,:) = [0 0 zmax + leUnb];
            end
        end
    end
    dofsGl = [dofsGl;dofsUn(1:iNew,:)];                                     % append to global dofs
    conne  = [conne;connUnb(1:iNew,:)];                                     % append to connectivity
    matNo  = [matNo;matNoUnb];                                              % append to material number
    sizes  = [sizes;ones(iNew,1)*leUnb];                                    % append to sizes
    coordV = [coordV;coordVunb];                                            % append to coordinates
    unbBTa = zeros(size(sizes));                                            % extend flags to include all elements                                         
    unbBTa(end-iNew+1 : end) = unbBT(1:iNew);                               
    isUnbounded = false(size(sizes));                                       % keep track of which elements are unbounded
    isUnbounded(end-iNew+1 : end) = true;
    eleNosUnb = nEle + (1:iNew)';
    for j = 1:nPde                                                          % loop PDEs
        pde2ele{j} = [pde2ele{j};eleNosUnb(unbPde==j)];                     % add new element numbers
    end
else                                                                        % no exact boundary conditions
    isUnbounded = false(size(sizes));                                       % no unbounded domains
    unbBTa = false;
end


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
