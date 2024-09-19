%% check basic properties of coupling
function [typeC,cOrder,kBotInd,kTopInd,nParameters,R] = analyzeCoupling(obj,R)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
Rp1 = R{1};
Rp2 = R{2};
Rs1 = R{3};
Rs2 = R{4};

waveVel = obj.waveVel;                                                      % wave velocities in unbounded media

existMatrices = ~[cellfun(@isempty,obj.Rp);cellfun(@isempty,obj.Rs)]';      % which matrices exist: [Rp,bot Rs,bot; Rp,top Rs,top]
if sum(existMatrices(:)) == 4                                               % all matrices present
    typeC = 'SS';                                                           % two solid unbounded domains
    cOrder = waveVel([1 3 2 4]);                                            % wave velocities in the order needed in solver
    R = {Rp1,Rs1,Rp2,Rs2};                                                  % coupling matrices to pass to solver
    kBotInd = [1,2];                                                        % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
    kTopInd = [3,4];                                                        % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    nParameters = 6;                                                        % number of parameters in multiparameter eigenvalue problem
elseif sum(existMatrices(:)) == 3                                           % three matrices present
    typeC = 'FS';                                                           % one fluid, one solid
    if sum(existMatrices(1,:))==2                                           % first is solid
        cOrder = waveVel([1 3 2]);                                          % wave velocities in the order needed in solver
        R = {Rp1,Rs1,Rp2};                                                  % coupling matrices to pass to solver
        kBotInd = [1,2];                                                    % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
        kTopInd = 3;                                                        % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    else
        cOrder = waveVel([2,4,1]);                                          % wave velocities in the order needed in solver
        R = {Rp2,Rs2,Rp1};                                                  % coupling matrices to pass to solver
        kBotInd = 3;                                                        % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
        kTopInd = [1,2];                                                    % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    end
    nParameters = 5;                                                        % number of parameters in multiparameter eigenvalue problem
elseif sum(existMatrices(:)) == 1                                           % one matrix present
    typeC = 'F';                                                            % one fluid
    if existMatrices(1,1)                                                   % bottom
        cOrder = waveVel(1);                                                % wave velocities in the order needed in solver
        R = {Rp1};
        kBotInd = 1;                                                        % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
        kTopInd = [];                                                       % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    else                                                                    % top
        cOrder = waveVel(2);                                                % wave velocities in the order needed in solver
        R = {Rp2};
        kBotInd = [];                                                       % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
        kTopInd = 1;                                                        % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    end
    nParameters = 3;                                                        % number of parameters in multiparameter eigenvalue problem
elseif sum(existMatrices(:,2)) == 1                                         % two matrices, one of them shear waves
    typeC = 'S';                                                            % one solid
    if existMatrices(1,1)                                                   % bottom
        cOrder = waveVel(1,:);                                              % wave velocities in the order needed in solver
        R = {Rp1,Rs1};
        kBotInd = [1,2];                                                    % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
        kTopInd = [];                                                       % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    else                                                                    % top
        cOrder = waveVel(2,:);                                              % wave velocities in the order needed in solver
        R = {Rp2,Rs2};
        kBotInd = [];                                                       % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
        kTopInd = [1,2];                                                    % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    end
    nParameters = 4;                                                        % number of parameters in multiparameter eigenvalue problem
elseif sum(existMatrices(:,2)) == 0                                         % two matrices, none of them shear waves
    typeC = 'FF';                                                           % two fluids
    cOrder = waveVel(:,1).';                                                  % wave velocities in the order needed in solver
    R = {Rp1,Rp2};
    kBotInd = 1;                                                            % indices of eigenvalue solutions corresponding to wavenumbers at lower surface
    kTopInd = 2;                                                            % indices of eigenvalue solutions corresponding to wavenumbers at upper surface
    nParameters = 4;                                                        % number of parameters in multiparameter eigenvalue problem
end

% NOTE: If both sides are coupled to the same material, the coupling terms
% are all included in Rp1 and (for solid) Rs1


end


