

function obj = loopFrequencies(obj,nParameters,dofs,computeEigenvectors,kBotInd,kTopInd,cOrder,coeffM,R,typeC,dofD,aDofs,~,wbName)

% read matrices
% For now, we use full matrices in this solver, as it only works for plates
% and cylinders anyway
E0  = full(coeffM{1});
E11 = full(coeffM{2});
E12 = full(coeffM{3});
E2  = full(coeffM{4});
M0  = full(coeffM{6});

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% translate into common notation, adjust later
L1 =  E12 - E11;
L2 =  E0;
L0 = -E2;

nSols = 2^nParameters*dofs;                                                 % (maximum) number of possible solutions

f = linspace(obj.fMin,obj.fMax,obj.nSteps);                                 % frequency array
w = f*2*pi;                                                                 % circular frequencies
nF = length(f);                                                             % number of frequencies

% allocate results
k   = nan(nSols,nF);                                                        % wavenumbers
att = nan(nSols,nF);                                                        % attenuation
cp  = nan(nSols,nF);                                                        % phase velocities
if computeEigenvectors
    cg  = nan(nSols,nF);                                                    % group velocities
    phi = zeros(dofs,nSols,nF);                                             % eigenvectors of linearized problem
else
    cg  = [];
    phi = [];
end
kz_bottom = zeros(nSols,nF,numel(kBotInd));                                 % allocate wavenumbers in unbounded domain
kz_top    = zeros(nSols,nF,numel(kTopInd));

maxKfound = 0;                                                              % keep track of maximum number of solutions found
disp('computing dispersion curves...')
timeStart = tic;
for i = 1:nF

    kF = w(i)./cOrder;                                                      % free-field wave numbers in unbounded medium

    if computeEigenvectors
        [ki, kAll, vi] = eig_Leaky_all(L2,L1,L0,M0,R,typeC,kF,w(i),dofD);   % call general multipareig solver
        vi = [vi{:,1}];
    else
        [ki, kAll    ] = eig_Leaky_all(L2,L1,L0,M0,R,typeC,kF,w(i),dofD);   % call general multipareig solver
    end
    nki = numel(ki);                                                        % number of eigenvalues

    [ki, ind] = sort(ki);
    kAll = kAll(ind,:);

    kz_bottom(1:nki,i,:) = kAll(:,kBotInd+1);                               % store wavenumbers in unbounded domains
    kz_top(1:nki,i,:)    = kAll(:,kTopInd+1);
    k(1:nki,i) = ki;                                                        % store wavenumber
    att(1:nki,i) = imag(ki)*20/log(10)*1000;                                % store attenuation
    if computeEigenvectors                                                  % store eigenvectors
        vi = vi(:,ind);
        phi(aDofs,1:nki,i) = vi;
    end
    cp(1:nki,i) = w(i)./real(ki);                                           % store phase velocity

    maxKfound = max(maxKfound,nki);

    if (floor(mod(i,nF/100)) == 0) || (i==nF)
        wbName = myWaitbar(i,nF,wbName,'value');
    end

end


if obj.sameUnbMedia
    kz_top = -kz_bottom;
end

kz_top    = squeeze(kz_top);
kz_bottom = squeeze(kz_bottom);

timeEnd = toc(timeStart);
disp(['   ... finished in ',num2str(round(timeEnd,2)), ' s'])

% store with solver object
obj.f     = f;
obj.omega = w;
obj.k     = k(1:maxKfound,:);
obj.att   = att(1:maxKfound,:);
obj.cp    = cp(1:maxKfound,:);
if computeEigenvectors
    obj.cg = cg(1:maxKfound,:);
    obj.phi = phi(:,1:maxKfound,:);
    obj.u = obj.phi;
end
% obj.q = phi((1:end/2)+end/2,:,:);
obj.kz_bottom = kz_bottom(1:maxKfound,:,:);
obj.kz_top    = kz_top(1:maxKfound,:,:);
obj.tCPU      = timeEnd;

end
