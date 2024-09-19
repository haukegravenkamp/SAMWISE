function [obj,bcd,glb] = computeSolution(obj,glb,msh,geo,bcd,mat,opt)

%% interface
computeEigenvectors       = opt.numerics.computeEigenvectors;
unboundedRemoveIncoming   = opt.numerics.unboundedRemoveIncoming;
removeNegativeAttenuation = opt.numerics.removeNegativeAttenuation;
matLayers = mat(msh.materialNumber);

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% waitbar
wbName = initializeWaitbar(1);

%% create coupling matrices describing unbounded domains
[obj, glb, bcd] = getCouplingMatrices(obj,msh,glb,geo,bcd,mat);

Rp1 = obj.Rp{1};                                                            % bottom, pressure waves
Rp2 = obj.Rp{2};                                                            % top,    pressure waves
Rs1 = obj.Rs{1};                                                            % bottom, shear waves
Rs2 = obj.Rs{2};                                                            % top,    shear waves

RO = {Rp1,Rp2,Rs1,Rs2};                                                     % collect in one cell

%% get coefficient matrices

% quadratic eigenvalue problem without boundary conditions:
% [-k²*E0 + ik*(E12-E11) - E2 + w*C0 + w²*M0 ]*u = 0
coeffMO = readCoeffMatrices(glb);

%% fix dofs according to boundary conditions
dofs = size(coeffMO{1},1);
fixedDofs = [glb.fixedDofs,false(1,dofs-numel(glb.fixedDofs))];             % fixed dofs
aDofs = ~ fixedDofs;                                                        % active dofs
[coeffM,R] = selectMatrixIndices(aDofs,coeffMO,RO);                         % remove fixed dofs

%% check basic properties of coupling
[typeC,cOrder,kBotInd,kTopInd,nParameters,R] = analyzeCoupling(obj,R);

%% check if problem can be efficiently block-decomposed
[coeffM,R,dofD] = blockDecomposition(obj,coeffM,R,typeC,opt);

%% loop frequencies
obj = loopFrequencies(obj,nParameters,dofs,computeEigenvectors,kBotInd,kTopInd,cOrder,coeffM,R,typeC,dofD,aDofs,matLayers,wbName);

%% compute vertical power flux at surfaces
obj = computePowerFlux(obj,bcd,mat);

%% filter solutions
if removeNegativeAttenuation
    obj = removeNegativeAtt(obj);
end

if unboundedRemoveIncoming
    obj = removeIncomingWaves(obj,opt);
end

obj = sortEigenvalues(obj);

%% check residual
if opt.numerics.computeResidual
    E0  = coeffM{1};
    E11 = coeffM{2};
    E12 = coeffM{3};
    E2  = coeffM{4};
    M0  = coeffM{6};
    obj = checkResidual(obj,E0,E12-E11,E2,M0,typeC,R,aDofs);
end


end
