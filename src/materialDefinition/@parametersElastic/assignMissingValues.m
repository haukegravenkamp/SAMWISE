
function obj = assignMissingValues(obj,~)

%% interface
isotropy = obj.isotropy;
if all(size(obj.D)==[6 6]) && ~isempty(obj.rho)                             % elasticity matrix and density already set
    isotropy = -1;                                                          % code for "do nothing"
end
% fill elastic parameters depending on which are given
switch isotropy                                                             % directional dependency (isotropic, orthotropic etc)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
    case -1                                                                 % unknown but D matrix already set

        if isempty(obj.cs)                                                  % shear wave velocity missing (may be required for automatic element order)
            obj.cs = sqrt(min(eig(obj.D))/obj.rho);                         % compute based on smallest eigenvalue
        end

    case 0                                                                  % isotropic

        if isSet(obj,{'E','rho','nu'})
            [obj, isHandle] = makeHandles(obj,{'E','rho','nu'});
            if isHandle
                obj.G      = @(y,z) obj.E(y,z)./(2*(1+obj.nu(y,z)));
                obj.cs     = @(y,z) sqrt(obj.G(y,z)./obj.rho(y,z));
                obj.cl     = @(y,z) sqrt(2*obj.G(y,z)./obj.rho(y,z).*(1-obj.nu(y,z))./(1-2*obj.nu(y,z)));
                obj.lambda = @(y,z) 2*obj.G(y,z).*obj.nu(y,z)./(1-2*obj.nu(y,z));
            else
                obj.G      = obj.E./(2*(1+obj.nu));
                obj.cs     = sqrt(obj.G./obj.rho);
                obj.cl     = sqrt(2*obj.G./obj.rho.*(1-obj.nu)./(1-2*obj.nu));
                obj.lambda = 2*obj.G.*obj.nu./(1-2*obj.nu);
            end
        elseif isSet(obj,{'G','rho','nu'})
            [obj, isHandle] = makeHandles(obj,{'G','rho','nu'});
            if isHandle
                obj.cs     = @(y,z) sqrt(obj.G(y,z)./obj.rho(y,z));
                obj.cl     = @(y,z) sqrt(2*obj.G(y,z)./obj.rho(y,z).*(1-obj.nu(y,z))./(1-2*obj.nu(y,z)));
                obj.E      = @(y,z) 2*obj.G(y,z).*(1+obj.nu(y,z));
                obj.lambda = @(y,z) 2*obj.G(y,z).*obj.nu(y,z)./(1-2*obj.nu(y,z));
            else
                obj.cs     = sqrt(obj.G./obj.rho);
                obj.cl     = sqrt(2*obj.G./obj.rho.*(1-obj.nu)./(1-2*obj.nu));
                obj.E      = 2*obj.G.*(1+obj.nu);
                obj.lambda = 2*obj.G.*obj.nu./(1-2*obj.nu);
            end

        elseif isSet(obj,{'E','G','rho'})
            [obj, isHandle] = makeHandles(obj,{'E','G','rho'});
            if isHandle
                obj.nu     = @(y,z) obj.E(y,z)/2/obj.G(y,z)-1;
                obj.cs     = @(y,z) sqrt(obj.G(y,z)./obj.rho(y,z));
                obj.cl     = @(y,z) sqrt(2*obj.G(y,z)./obj.rho(y,z).*(1-obj.nu(y,z))./(1-2*obj.nu(y,z)));
                obj.lambda = @(y,z) 2*obj.G(y,z).*obj.nu(y,z)./(1-2*obj.nu(y,z));
            else
                obj.nu     = obj.E/2/obj.G-1;
                obj.cs     = sqrt(obj.G./obj.rho);
                obj.cl     = sqrt(2*obj.G./obj.rho.*(1-obj.nu)./(1-2*obj.nu));
                obj.lambda = 2*obj.G.*obj.nu./(1-2*obj.nu);
            end

        elseif isSet(obj,{'G','lambda','rho'})
            [obj, isHandle] = makeHandles(obj,{'G','lambda','rho'});
            if isHandle
                obj.nu = @(y,z) obj.lambda(y,z)/2/(obj.lambda(y,z)+obj.G(y,z));
                obj.cs = @(y,z) sqrt(obj.G(y,z)/obj.rho(y,z));
                obj.cl = @(y,z) sqrt(2*obj.G(y,z)/obj.rho(y,z)*(1-obj.nu(y,z))/(1-2*obj.nu(y,z)));
                obj.E  = @(y,z) 2*obj.G(y,z)*(1+obj.nu(y,z));
            else
                obj.nu = obj.lambda/2/(obj.lambda+obj.G);
                obj.cs = sqrt(obj.G/obj.rho);
                obj.cl = sqrt(2*obj.G/obj.rho*(1-obj.nu)/(1-2*obj.nu));
                obj.E  = 2*obj.G*(1+obj.nu);
            end

        elseif isSet(obj,{'cs','cl','rho'})
            [obj, isHandle] = makeHandles(obj,{'cs','cl','rho'});
            if isHandle
                obj.nu     = @(y,z) (obj.cl(y,z)^2/2-obj.cs(y,z)^2)/(obj.cl(y,z)^2-obj.cs(y,z)^2);
                obj.G      = @(y,z) obj.cs(y,z)^2*obj.rho(y,z);
                obj.E      = @(y,z) 2*obj.G(y,z)*(1+obj.nu(y,z));
                obj.lambda = @(y,z) 2*obj.G(y,z)*obj.nu(y,z)/(1-2*obj.nu(y,z));
            else
                obj.nu     = (obj.cl^2/2-obj.cs^2)/(obj.cl^2-obj.cs^2);
                obj.G      = obj.cs^2*obj.rho;
                obj.E      = 2*obj.G*(1+obj.nu);
                obj.lambda = 2*obj.G*obj.nu/(1-2*obj.nu);
            end
        elseif isSet(obj,{'cs','nu','rho'})
            [obj, isHandle] = makeHandles(obj,{'cs','nu','rho'});
            if isHandle
                obj.G      = @(y,z) obj.cs(y,z)^2*obj.rho(y,z);
                obj.E      = @(y,z) 2*obj.G(y,z)*(1+obj.nu(y,z));
                obj.lambda = @(y,z) 2*obj.G(y,z)*obj.nu(y,z)/(1-2*obj.nu(y,z));
                obj.cl     = @(y,z) sqrt(2*obj.G(y,z)/obj.rho(y,z)*(1-obj.nu(y,z))/(1-2*obj.nu(y,z)));
            else
                obj.G      = obj.cs^2*obj.rho;
                obj.E      = 2*obj.G*(1+obj.nu);
                obj.lambda = 2*obj.G*obj.nu/(1-2*obj.nu);
                obj.cl     = sqrt(2*obj.G/obj.rho*(1-obj.nu)/(1-2*obj.nu));
            end
        elseif isSet(obj,{'cl','nu','rho'})
            [obj, isHandle] = makeHandles(obj,{'cl','nu','rho'});
            if isHandle
                obj.cs     = @(y,z) sqrt(obj.cl(y,z)^2*(1/2-obj.nu(y,z))/(1-obj.nu(y,z)));
                obj.G      = @(y,z) obj.cs(y,z)^2*obj.rho(y,z);
                obj.E      = @(y,z) 2*obj.G(y,z)*(1+obj.nu(y,z));
                obj.lambda = @(y,z) 2*obj.G(y,z)*obj.nu(y,z)/(1-2*obj.nu(y,z));
            else
                obj.cs     = sqrt(obj.cl^2*(1/2-obj.nu)/(1-obj.nu));
                obj.G      = obj.cs^2*obj.rho;
                obj.E      = 2*obj.G*(1+obj.nu);
                obj.lambda = 2*obj.G*obj.nu/(1-2*obj.nu);
            end
        else
            error('not enough elastic parameters provided')
        end

        obj.D = obj.computeD_isotropic(obj.G, obj.nu, 3, isHandle);         % compute isotropic elasticity matrix

    case 1                                                                  % transversely isotropic
        error('material behavior not implemented')
    case 2                                                                  % orthotropic
        % Works only if properties are aligned with coordinate system
        E1   = obj.E1;   E2   = obj.E2;   E3   = obj.E3;
        nu12 = obj.nu12; nu13 = obj.nu13; nu23 = obj.nu23;
        G23  = obj.G23;  G31  = obj.G31;  G12  = obj.G12;
        % D-matrix
        obj.D = parametersElastic.computeD_orthotropic(...
            E1,E2,E3, nu12,nu13,nu23, G23,G31,G12);

    case 3                                                                  % anisotropic
        % assumes the required elasticity matrix D is
        % provided directly -> do nothing
end


end
