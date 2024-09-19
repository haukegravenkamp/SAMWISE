function [obj,glb,bcd] = getCouplingMatrices(obj,msh,glb,geo,bcd,mat)


%% check if boundary conditions for unbounded medium exist
unbType = zeros(numel(bcd),1);                                              % which boundary conditions for unbounded media exist
% 0: none, 1: fluid coupling, 2: solid coupling
isUnbounded = false(numel(bcd),1);                                          % which boundary conditions are for unbounded media
botTop = zeros(numel(bcd),1);                                               % bottom or top
for ib = 1:numel(bcd)                                                       % loop boundary conditions

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
    if isa(bcd(ib),'bcUnbounded')                                           % is a b.c. for unbounded domain
        matNoUnb = bcd(ib).material;                                        % material number of unbounded medium
        isUnbounded(ib) = true;                                             % remember bc number
        if strcmp(mat(matNoUnb).behavior,'acoustic')                        % check whether coupling is to acoustic or elastic medium
            unbType(ib) = 1;
        elseif strcmp(mat(matNoUnb).behavior,'elastic')
            unbType(ib) = 2;
        end
    end
    if strcmp(bcd(ib).location,'bottom')
        botTop(ib) = 1;
    elseif strcmp(bcd(ib).location,'top')
        botTop(ib) = 2;
    end

end
if ~any(unbType)                                                            % no unbouded boundary conditions found
    error(['Solver for embedded plates called but no embedded boundary' ...
        ' conditions found. Check for missing boundary conditions ' ...
        'or use standard solver.'])
end

bcUnb = bcd(isUnbounded);                                                   % unbounded boundary conditions
unbInd = find(isUnbounded);                                                 % indices of unbounded bc in global sorting
unbType = unbType(isUnbounded);
botTop = botTop(isUnbounded);
% check number of dofs per node in elastic problem
if strcmp(geo.assumption2D,'planeStrain') || ...
        strcmp(geo.assumption2D,'planeStress')
    nDofsElastic = 2;
else
    nDofsElastic = 3;
end
nDofsBC = sum(unbType==1) + sum(unbType==2)*nDofsElastic;                   % total number of dofs to be added for b.c.
nDofsFree = msh.nDofs;                                                      % number of dofs without coupling
nDofsCoupled = nDofsFree + nDofsBC;                                         % updated number of dofs
dofsUnb = nan(2,3);                                                         % allocate dofs of unbounded domains
dofsUnb(unbType==1,1) = 1;                                                  % dofs of acoustic coupling
dofsUnb(unbType==2,1:nDofsElastic)=repmat(1:nDofsElastic,sum(unbType==2),1);% dofs of elastic coupling
dofsUnb(2,:) = dofsUnb(2,:) + nDofsFree + numel(nonnan(dofsUnb(1,:)));      % add previous dofs
dofsUnb(1,:) = dofsUnb(1,:) + nDofsFree;
dofsUnb = dofsUnb(isUnbounded,:);                                           % remove unneeded rows

% store dofs with original boundary condition array
bcInd = 1:numel(bcd);
bcInd = bcInd(isUnbounded);
for i = 1:numel(bcInd)
    bcd(bcInd(i)).dofsUnbounded = nonnan(dofsUnb(i,:));
end

%% read existing matrices and extend
E0  = glb.feMatrices{1,1};
E11 = glb.feMatrices{2,1};
E12 = glb.feMatrices{3,1};
E2  = glb.feMatrices{4,1};
C0  = glb.feMatrices{1,2};
M0  = glb.feMatrices{1,3};
E0(nDofsCoupled,nDofsCoupled)  = 0;
E11(nDofsCoupled,nDofsCoupled) = 0;
E12(nDofsCoupled,nDofsCoupled) = 0;
E2(nDofsCoupled,nDofsCoupled)  = 0;
C0(nDofsCoupled,nDofsCoupled)  = 0;
M0(nDofsCoupled,nDofsCoupled)  = 0;

%% allocate coupling matrices
Rp = {[],[]};                                                               % pressure waves, bottom/top
Rs = {[],[]};                                                               % shear waves, bottom/top
if any(botTop==1)                                                           % if coupling at bottom exists
    Rp{1} = zeros(nDofsCoupled);                                            % allocate bottom, pressure wave
    if unbType(botTop==1)==2                                                % if it's elastic coupling
        Rs{1} = zeros(nDofsCoupled);                                        % bottom, shear wave
    end
end
if any(botTop==2)                                                           % if coupling at top exists
    Rp{2} = zeros(nDofsCoupled);                                            % top, pressure wave
    if unbType(botTop==2)==2                                                % elastic coupling
        Rs{2} = zeros(nDofsCoupled);                                        % top, shear wave
    end
end


%% compute coupling terms
waveVel = zeros(2);                                                         % for storing wave velocities in unbounded medium

for ib = 1:numel(bcUnb)                                                     % loop unbouded boundary conditions

    dofsEle = bcUnb(ib).dofs;                                               % dofs corresponding to boundary condition in element
    dofsUnbC = nonnan(dofsUnb(ib,:));
    matNoUnb = bcUnb(ib).material;                                          % material number of unbounded medium
    cl  = mat(matNoUnb).parameters.cl;                                      % longitudinal wave velocity

    if strcmp(bcUnb(ib).location,'bottom')                                  % bottom of plate structure
        signBC = -1;                                                        % invert sign of tractions
        matInd = botTop(ib);                                                % index of coupling matrices
    else
        signBC = 1;
        matInd = botTop(ib);
    end
    rho = mat(matNoUnb).parameters.rho;                                     % density

    if unbType(ib) == 1                                                     % acoustic coupling

        dofsEle = dofsEle(end);                                               % vertical dof

        % bottom: pressure causes positive traction, negative on lhs
        % top:    pressure causes negative traction, positive on lhs
        E2(dofsEle,dofsUnbC) = signBC;                                      % condition on solid (bottom): n*sigma(y=0) = p         -

        % bottom: positive acceleration causes negative pressure gradient, positive
        % in frequency domain, negative on lhs
        % top: positive acceleration causes positive pressure gradient, negative
        % in frequency domain, positive on lhs
        M0(dofsUnbC,dofsEle) = signBC*rho;                                  % pressure equation (bottom): dp/dn = -rhoF * n*u'' =  omega^2*rhoF*u   -
        Rp{matInd}(dofsUnbC,dofsUnbC) = -signBC;

        cs = 0;

    elseif unbType(ib) == 2                                                 % elastic coupling

        lbd = mat(matNoUnb).parameters.lambda;                              % Lamé parameter lambda
        mu  = mat(matNoUnb).parameters.G;                                   % Lamé parameter mu
        
        cs  = mat(matNoUnb).parameters.cs;                                  % shear wave velocity
        % Dxy = [0 mu 0; lbd 0 0; 0 0 0];
        % Dyy = [mu 0 0; 0 lbd+2*mu 0; 0 0 mu];

        Dxz = [0 0 mu; 0 0 0; lbd 0 0];
        Dzz = [mu 0 0; 0 mu 0; 0 0 lbd+2*mu];

        % helper matrices
        D1 = diag([1 0 0]);
        D2 = diag([0 1 1]);
        % A0 = [1 0 0; 0 -1 0; 0 0 1];
        % A1 = [0 0 0; 1 0 0; 0 0 0];
        % A2 = [0 1 0; 0 0 0; 0 0 0];
        A0 = [1 0 0; 0 1 0; 0 0 -1];
        A1 = [0 0 0; 0 0 0; 1 0 0];
        A2 = [0 0 1; 0 0 0; 0 0 0];
        nDofsEle = numel(dofsEle);
        if nDofsEle==2                                                % plane strain
            ind = [1,3];
            Dxz = Dxz(ind,ind);
            Dzz = Dzz(ind,ind);
            D1 = D1(ind,ind);
            D2 = D2(ind,ind);
            A0 = A0(ind,ind);
            A1 = A1(ind,ind);
            A2 = A2(ind,ind);
        end

        % signBC = -signBC;

        r0 = -(Dxz*A0 - Dzz*(A1*D1+A2*D2));
        r1 = -(Dxz*A1 + Dzz*A0*D1);
        r2 = -(Dxz*A2 + Dzz*A0*D2);
        r3 = -Dzz*A1*D1;
        r4 = -Dzz*A2*D2;

        % contributions to equation for unbounded medium 
        E11(dofsUnbC,dofsEle)         =  eye(nDofsEle);
        E0(dofsUnbC,dofsUnbC)         =  A0;
        Rp{matInd}(dofsUnbC,dofsUnbC) = -A1;
        Rs{matInd}(dofsUnbC,dofsUnbC) = -A2;

        % contributions to plate equation
        E0(dofsEle,dofsUnbC)         = -r0*signBC;
        Rp{matInd}(dofsEle,dofsUnbC) =  r1*signBC;
        Rs{matInd}(dofsEle,dofsUnbC) =  r2*signBC;
        M0(dofsEle,dofsUnbC)         =  r3/cl^2*signBC + r4/cs^2*signBC;
        % paper version: assembled into E2 instead of M0 -> opposite sign 

        bcd(unbInd(ib)).rMatrices = {r0,r1,r2,r3,r4};
        
    end

    waveVel(matInd,:)=[cl,cs];                                              % store wave velocities

end

if numel(bcUnb)==2 && (bcUnb(1).material == bcUnb(2).material)              % b.c. at top and bottom and same material
    Rp{1} = Rp{1} - Rp{2};                                                  % add coupling matrices and remove second one
    Rp{2} = [];                                                             % signs are treated as if unbounded is coupled to bottom side
    Rs{1} = Rs{1} - Rs{2};                                                  % i.e., negative wavenumbers propagate away from plate
    Rs{2} = [];
    sameUnbMedia = true;
else
    sameUnbMedia = false;
end


%% output
glb.feMatrices{1,1} = E0;
glb.feMatrices{2,1} = E11;
glb.feMatrices{3,1} = E12;
glb.feMatrices{4,1} = E2;
glb.feMatrices{1,2} = C0;
glb.feMatrices{1,3} = M0;
obj.Rp = Rp;
obj.Rs = Rs;
obj.waveVel = waveVel;
obj.sameUnbMedia = sameUnbMedia;


end

