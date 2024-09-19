function glb = dashpotMatrix(obj,glb,bcd,mat)


%% interface
C0 = glb.feMatrices{1,2};                                                   % damping matrix, linear in 1i*omega, no spatial derivative

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% add terms to damping matrix for each boundary condition
% if there are several b.c. the terms would be added up

isCylinder = isa(obj,'mesh1Dcylinder');                                     % whether geometry is a cylinder

if isempty(C0)
    C0 = zeros(obj.nDofs);
end
for ib = 1:numel(bcd)                                                       % loop boundary conditions
    
    dofs = bcd(ib).dofs;                                                    % dofs corresponding to boundary condition
    matNoUnb = bcd(ib).material;                                            % material number of unbounded medium
    matNoEle = bcd.materialNo;                                              % material number of element 
    
    if isCylinder
        r0 = bcd.coord(end);
    else
        r0 = 1;
    end

    cl = mat(matNoUnb).parameters.cl;                                       % longitudinal wave velocity of unbounded medium
    if isprop(mat(matNoUnb).parameters,'cs')                                % if unbounded medium has a shear velocity
        cs = mat(matNoUnb).parameters.cs;                                   % read shear velocity
        isElasticUnb = true;                                                % assume this material to be elastic
    else
        isElasticUnb = false;                                               % assume this material to be acoustic
    end
    rho = mat(matNoUnb).parameters.rho;                                     % mass density of unbounded medium
    if isprop(mat(matNoEle).parameters,'cs')                                % if element material has a shear velocity
        isElasticEle = true;                                                % assume this material to be elastic
    else
        isElasticEle = false;                                               % assume this material to be acoustic
    end

    nDofs = numel(dofs);                                                    % number of dofs affected by boundary condition

    if isElasticEle && isElasticUnb                                         % elastic element, elastic unbounded
        if nDofs == 2
            cAll = [cs,cl];
        else
            cAll = [cs,cs,cl];                                              % three dofs coupling matrix
        end
        C0(dofs,dofs) = C0(dofs,dofs) - rho*r0*diag(cAll(1:nDofs));         % use adequate entries in coupling matrix
    elseif isElasticEle && ~isElasticUnb                                    % elastic element, acoustic unbounded
        C0(dofs(end),dofs(end)) = C0(dofs(end),dofs(end)) -rho*r0*cl;       % couple only longitudinal wave to vertical component
    else                                                                    % element is fluid
        C0(dofs,dofs) = C0(dofs,dofs) -rho*r0*cl;                           % couple only longitudinal wave
    end

end

%% output
glb.feMatrices{1,2} = C0;

end
