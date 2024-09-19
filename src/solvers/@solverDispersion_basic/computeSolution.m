function [obj, bcd, glb] = computeSolution(obj,glb,msh,geo,bcd,mat,opt)
                                     
disp('computing dispersion curves...')
timeStart = tic;

if (msh.nDofs > 200)                                                       % expensive problems
    showWaitbars = true;                                                    % show waitbars while solving
else                                                                        % small problems
    showWaitbars = [];                                                      % don't bother showing progress
end
wbName = initializeWaitbar(showWaitbars);


% the basic solver assumes the existence of the standard matrices as known
% for acoustic and elastic media. Can be generalized if needed.
% The quadratic eigenvalue problem reads
% [-k²*E0 + ik*(E12-E11) - E2 + w*C0 + w²*M0 ]*u = 0
E0  = glb.feMatrices{1,1};
E11 = glb.feMatrices{2,1};
E12 = glb.feMatrices{3,1};
E2  = glb.feMatrices{4,1};
C0  = glb.feMatrices{1,2};
if isempty(C0)
    C0 = 0*E0;
end
M0  = glb.feMatrices{1,3};
if isempty(M0)
    M0 = 0*E0;
end
dofs = size(E0,1);

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

nSolutions = obj.nSolutions;                                                % number of solutions requested by user
nEigenvaluesDefault = opt.numerics.nEigenvaluesDefault;                              % default number of solutions
useSparseSize = opt.numerics.useSparseSize;

% set number of solutions to compute
if isscalar(nSolutions)                                                     % number of solutions given explicitly
    nSolutions = abs(round(nSolutions));                                    % make sure it's a positive integer
elseif ~isscalar(nSolutions)&&strcmp(nSolutions,'all')                      % all solutions requested
    nSolutions = dofs;                                                      % compute all solutions
else                                                                        % automatic setting
    if ~ischar(nSolutions)||~strcmp(nSolutions,'auto')                      % it was none of the known inputs
        warning('Unknown value provided as solver setting "nSolutions". Using defaults.')
    end
    if dofs <= useSparseSize                                                % small matrices
        nSolutions = dofs;                                                  % compute all solutions
    else                                                                    % large matrices
        nSolutions = nEigenvaluesDefault;                                   % default
    end
end

if dofs<=useSparseSize
    E0 = full(E0);
    E11 = full(E11);
    E12 = full(E12);
    E2 = full(E2);
    M0 = full(M0);
    C0 = full(C0);
    sparseMatrices = false;
else
    sparseMatrices = true;
end

fixedDofs = glb.fixedDofs;                                                  % fixed dofs
aDofs = ~fixedDofs;                                                        % active dofs
indPhi = [aDofs,aDofs].';
[E0,E11,E12,E2,C0,M0] = selectMatrixIndices(aDofs,E0,E11,E12,E2,C0,M0);   % remove fixed dofs


%% check if problem can be efficiently block-decomposed
decomposeEVP = opt.numerics.decomposeEVP;

if decomposeEVP
    [ind,nBl] = obj.blockDecomposition(E2+M0,1e-6);
    E1 = E12-E11;
    AC = E0 + E2 + M0 + C0;
end
decomposed = false;
if decomposeEVP && (nBl == 2) ...                                           % two blocks (plane strain or plane stress)
        && ~any(E1([ind{1}],[ind{1}]),'all') ...
        && ~any(E1([ind{2}],[ind{2}]),'all') ...
        && ~any(AC([ind{1}],[ind{2}]),'all')...
        && ~any(AC([ind{2}],[ind{1}]),'all')

    i1 = ind{1};
    i2 = ind{2};
    iA = [i1,i2];

    Md = M0(iA,iA);
    Cd = C0(iA,iA);

    AB = [ E0(i1,iA); E1(i2,i1), E0(i2,i2)];
    BC = [-E2(i1,i1), E1(i1,i2);-E2(i2,iA)];

    [~,dofD]=sort(iA);                                                      % indices to revert decomposition

    decomposed = true;

elseif decomposeEVP && (nBl == 3) ...                                       % three blocks (including out-of-plane modes)
        && ~any(E1([ind{1}],[ind{1}]),'all') ...
        && ~any(E1([ind{2}],[ind{2}]),'all') ...
        && ~any(E1([ind{3}],[ind{3}]),'all') ...
        && ~any(AC([ind{1}],[ind{2}]),'all')...
        && ~any(AC([ind{1}],[ind{3}]),'all')...
        && ~any(AC([ind{2}],[ind{1}]),'all')...
        && ~any(AC([ind{2}],[ind{3}]),'all')...
        && ~any(AC([ind{3}],[ind{1}]),'all')...
        && ~any(AC([ind{3}],[ind{2}]),'all')

    if ~any(E1([ind{1}],[ind{:}]),'all')                                    % find the isolated block of out-of-plane modes and move it to the end for consistency
        i1 = ind{2};
        i2 = ind{3};
        i3 = ind{1};
    elseif ~any(E1([ind{2}],[ind{:}]),'all')
        i1 = ind{1};
        i2 = ind{3};
        i3 = ind{2};
    elseif ~any(E1([ind{3}],[ind{:}]),'all')
        i1 = ind{1};
        i2 = ind{2};
        i3 = ind{3};
    end

    iA = [i1,i2,i3];

    Md = M0(iA,iA);
    Cd = C0(iA,iA);

    AB = [ E0(i1,iA); E1(i2,i1), E0(i2,[i2,i3]); E0(i3,iA)];
    BC = [-E2(i1,i1), E1(i1,[i2,i3]);-E2(i2,iA);-E2(i3,iA)];

    [~,dofD]=sort(iA);                                                      % indices to revert decomposition

    decomposed = true;

end

if decomposed
    
    % [BC,indb,Md,Cd,AB]=mesh.reduceBandwidth(BC,Md,Cd,AB);

    iAB = AB\eye(size(AB));
    iABBC = -iAB*BC;
    iABM  = -iAB*Md;
    iABC  = iAB*Cd*1i;
else
    Zs = solver.Zmatrix(0,E0,E11,E12,E2,M0,C0);                             % static Z-matrix
    sizeZ = size(Zs,1);
    M = zeros(sizeZ);                                                       % expand M0 to size of Z for convenience
    M(end/2+1:end,1:end/2) = M0;
    C = zeros(sizeZ);                                                       % same for C0
    C(end/2+1:end,1:end/2) = C0;
end

%% loop frequencies

f = linspace(obj.fMin,obj.fMax,obj.nSteps);                                 % frequency array
omega = 2*pi*f;                                                             % circular frequencies
nF = length(f);

% allocate results
k   = nan(2*dofs,nF);                                                       % wavenumbers
att = nan(2*dofs,nF);                                                       % attenuation
cp  = nan(2*dofs,nF);                                                       % phase velocities
cg  = nan(2*dofs,nF);                                                       % group velocities
phi = zeros(2*dofs,nSolutions,nF);                                              % eigenvectors of linearized problem

for i = 1:nF
    if decomposed                                                           % linearization by decomposition, see Kausel                                     
        if sparseMatrices
            [vi,lam2] = eigs(BC + Md*omega(i)^2 - Cd*1i*omega(i),-AB,nSolutions,'smallestabs');              % solve eigenvalue problem
            lam2 = diag(lam2);
        else
            ABC = iABBC + iABM*omega(i)^2 + iABC*omega(i);
            [vi,lam2] = eig(ABC,'vector');                                  % solve eigenvalue problem
        end
        lam = sqrt(lam2);
        lambda = [lam;-lam];
        vu1 = vi(dofD,:);
        vu2 = vu1;
        vu1(i2,:) =  vu1(i2,:)./lam.';
        vu2(i2,:) = -vu2(i2,:)./lam.';
        vq1 =  E0*vu1*diag(lam) + E12*vu1;
        vq2 = -E0*vu2*diag(lam) + E12*vu2;
        vi = [vu1,vu2;vq1,vq2];

    else                                                                    % linearization using 'Z-matrix'                              
        Z = Zs + M*omega(i)^2 - C*1i*omega(i);                              % dynamic Z-matrix
        
        if sparseMatrices
            [vi,lambda] = eigs(-Z,[],nSolutions,'smallestabs');                                      % solve eigenvalue problem
            lambda = diag(lambda);
        else
            [vi,lambda] = eig(-Z,'vector');                                     % solve eigenvalue problem
        end

    end
    ki = lambda/1i;                                                         % wave numbers
    vu = vi(1:end/2,:);                                                     % displacements
    vq = vi(end/2+1:end,:);                                                 % forces
    
    vgs = obj.computeGroupVelocity(ki,vu,vq,M0,omega(i));

    % store results at current frequency
    nSols = numel(ki);                                                      % number of found solutions
    k(1:nSols,i)     = ki;
    att(1:nSols,i)   = imag(ki)*20/log(10)*1000;
    phi(indPhi,1:nSols,i) = vi;
    cg(1:nSols,i)    = real(vgs);
    cp(1:nSols,i)    = omega(i)./real(ki);

    if floor(mod(i,nF/100)) == 0
        wbName = myWaitbar(i,nF,wbName,'value');
    end

end

disp(['   ... finished in ',num2str(round(toc(timeStart),2)), ' s'])

% store with solver object
obj.f = f;
obj.omega = omega;
obj.k = k;
obj.att = att;
obj.cp = cp;
obj.cg = cg;
obj.phi = phi;
obj.u = phi(1:end/2,:,:);
obj.q = phi((1:end/2)+end/2,:,:);

end
