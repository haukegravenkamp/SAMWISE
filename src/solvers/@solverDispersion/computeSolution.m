function [sol, bcd, glb] = computeSolution(obj,glb,msh,geo,bcd,mat,opt)

%% interface
unboundedModel = opt.model.unboundedModel;

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% check if boundary conditions for unbounded medium exist
isUnbounded = false(numel(bcd),1);                                          % which boundary conditions are for unbounded media
for ib = 1:numel(bcd)                                                       % loop boundary conditions
    if isa(bcd(ib),'bcUnbounded')                                           % is a b.c. for unbounded domain
        isUnbounded(ib) = true;
    end
end
if any(isUnbounded)
    existUnbounded = true;                                                  % set flag
else
    existUnbounded = false;                                                  % set flag
end
%% check which solver is to be used for unbounded media
useExact = false;
useApprox = false;
if strcmp(unboundedModel,'exact')                                           % use exact boundary conditions
    useExact = true;
    if ~isa(geo,'plate')                                                    % not a plate
        useExact = false;                                                   % cannot use exact b.c.
        if existUnbounded                                                   % unbounded medium present
            warning(['Exact embedded boundary conditions are only ' ...     
                'implemented for plates. Switching to dashpot boundary.'])  % warn
        end
    end
elseif strcmp(unboundedModel,'approx')                                      % approximate sqrt terms in boundary conditions
    useApprox = true;
else                                                                        % use dashpot approximation if any
    if existUnbounded && ~strcmp(unboundedModel,'dashpot')                  % warn about unknown setting
        warning(['unknown option in model.unboundedModel, ' ...
            'using dashpot approximation'])
    end
end

%% create adequate solver
% can be extended if additional solvers for other situations are introduced
if existUnbounded && useExact                                               % unbounded domain exists and exact solver to be used
    sol = solverDispersion_embeddedExact;                                   % create special solver for embedded waveguides
elseif existUnbounded && useApprox                                          % unbounded domain exists and approx solver to be used
    sol = solverDispersion_embeddedApproximate;                             % create special solver for embedded waveguides
else                                                                        % no unbounded domain or dashpot approximation
    if existUnbounded                                                       % unbounded medium exists
        glb = dashpotMatrix(msh,glb,bcd(isUnbounded),mat);                  % set boundary conditions
    end
    sol = solverDispersion_basic;                                           % create basic solver
end

%% copy common properties to specialized solver
sol.fMin = obj.fMin;
sol.fMax = obj.fMax;
sol.nSteps = obj.nSteps;
sol.nSolutions = obj.nSolutions;

%% call solver
[sol, bcd, glb] = computeSolution(sol,glb,msh,geo,bcd,mat,opt);

%% compute vertical derivatives
sol = getDerivatives(sol,msh,glb);



end
