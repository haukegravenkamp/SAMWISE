

function obj = loopFrequencies(obj,nParameters,dofs,computeEigenvectors,kBotInd,kTopInd,cOrder,coeffM,R,typeC,dofD,aDofs,matLayers)

cAll = nan(numel(matLayers),2);
for i = 1:numel(matLayers)
    cAll(i,1) = matLayers(i).parameters.cl;
    cs = matLayers(i).parameters.cs;
    if ~isempty(cs)
        cAll(i,2) = cs;
    end
end

cAll = [cAll(:);cOrder(:)];

E0  = coeffM{1};
E11 = coeffM{2};
E12 = coeffM{3};
E2  = coeffM{4};
M0  = coeffM{6};
E1 =  E12 - E11;

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
ky_bottom = zeros(nSols,nF,numel(kBotInd));                                 % allocate wavenumbers in unbounded domain
ky_top    = zeros(nSols,nF,numel(kTopInd));

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
maxKfound = 0;                                                              % keep track of maximum number of solutions found
disp('computing dispersion curves...')
timeStart = tic;
for i = 1:nF

    kF = w(i)./cOrder;                                                      % free-field wave numbers in unbounded medium

    if computeEigenvectors
        [ki, kAll, vi] = eig_Leaky_approximate(E0,E1,E2,M0,R,typeC,kF,w(i),cAll);
        vi = [vi{:,1}];
    else
        [ki, kAll] = eig_Leaky_approximate(E0,E1,E2,M0,R,typeC,kF,w(i),cAll);
    end
    nki = numel(ki);                                                        % number of eigenvalues

    [ki, ind] = sort(ki);
    kAll = kAll(ind,:);

    ky_bottom(1:nki,i,:) = kAll(:,kBotInd+1);                               % store wavenumbers in unbounded domains
    ky_top(1:nki,i,:)    = kAll(:,kTopInd+1);
    k(1:nki,i) = ki;                                                        % store wavenumber
    att(1:nki,i) = imag(ki)*20/log(10)*1000;                                % store attenuation
    if computeEigenvectors                                                  % store eigenvectors
        vi = vi(:,ind);
        phi(aDofs,1:nki,i) = vi;
    end
    cp(1:nki,i) = w(i)./real(ki);                                           % store phase velocity

    maxKfound = max(maxKfound,nki);

end

if obj.sameUnbMedia
    ky_top = -ky_bottom;
end

ky_top    = squeeze(ky_top);
ky_bottom = squeeze(ky_bottom);

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
obj.ky_bottom = ky_bottom(1:maxKfound,:,:);
obj.ky_top = ky_top(1:maxKfound,:,:);
obj.tCPU = timeEnd;

end
