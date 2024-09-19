function [uxz,pxz,sigxz,Pxz] = pointValuesElements(uEle,wKmode,Jdet,x,eta,pSF,itg,dof,nx,nz,matEle,pdeName)

w  = wKmode(1);
kx = wKmode(2);

uxz = zeros(nz,nx,3);                                                       % allocate displacement
Pxz = zeros(nz,nx,3);                                                       % allocate Poynting vector

nU = numel(uEle);

if strcmp(pdeName,'pdeAcoustics')
    rho = matEle.parameters.rho;                                            % read mass density
    pxz = zeros(nz,nx,1);                                                   % allocate displacement
    sigxz = [];

    for iz = 1:nz                                                           % loop vertical position in element

        shp = shpLagrange(pSF+1,eta(iz),itg);                               % get shape functions and derivatives
        shp = [shp(:,1),shp(:,end),shp(:,2:end-1)];
        N   = kron(shp(1,:),eye(dof));                                      % shape function matrix for all dofs
        dN  = kron(shp(2,:),eye(dof));                                      % shape function derivatives for all dofs
        for ix = 1:nx                                                       % loop horizontal position

        end
        % only for acoustic problem: particle acceleration is related
        % to pressure:
        % a = -1/rho grad(p) -> u = 1/(omega²*rho) grad(p)
        for ix = 1:nx                                                       % loop horizontal position
            p = N*uEle * exp(1i*kx*x(ix));                                  % interpolate pressure
            pdx =  N*uEle * exp(1i*kx*x(ix))*1i*kx;                         % dp/dx
            pdz = dN*uEle * exp(1i*kx*x(ix)) /Jdet(iz);                     % dp/dy
            gradP = [pdx, pdz, 0];                                          % grad(p)
            u = 1/w^2/rho * gradP;                                          % displacement
            v = -1i*w*u;                                                    % velocity
            Pxz(iz,ix,:) = real(v'*p)/2;                                   % Poynting vector

            pxz(iz,ix,:) = p;
            uxz(iz,ix,:) = u;
        end
    end
end

if strcmp(pdeName,'pdeElasticity')
    D = matEle.elastic.D;                                                   % read elasticity matrix
    sigxz = zeros(nz,nx,6);                                                 % allocate stresses
    switch dof
        case 1                                                              % out of plane
            uEle = reshape([zeros(2,nU);uEle.'],[],1);                      % only third component is nonzero
        case 2                                                              % plane strain/stress
            unonzero = reshape(uEle,2,[]);
            uEle = reshape([unonzero(1,:); zeros(1,nU/2); unonzero(2,:)],[],1);       % components 1 and 3 nonzero
    end
    for iz = 1:nz                                                         % loop vertical position in element

        shp = shpLagrange(pSF+1,eta(iz),itg);                               % get shape functions and derivatives
        shp = [shp(:,1),shp(:,end),shp(:,2:end-1)];
        N   = kron(shp(1,:),eye(3));                                      % shape function matrix for all dofs
        dN   = kron(shp(2,:),eye(3));                                     % shape function derivatives for all dofs
        for ix = 1:nx                                                       % loop horizontal position
            u   = N*uEle * exp(1i*kx*x(ix));                                % interpolate solution, regardless of physical meaning
            udx =  N*uEle * exp(1i*kx*x(ix))*1i*kx;                         % du/dx
            udy = 0*uEle;                                                   % du/dy
            udz = dN*uEle * exp(1i*kx*x(ix)) /Jdet(iz);                     % du/dz
            e = [udx(1), udy(2), udz(3),...
                udy(3)+udz(2), udx(3)+udz(1), udy(1)+udx(2)].';             % strain
            sig = D*e;                                                      % stress (Voigt notation)
            sigT = voigt(sig);                                              % stress tensor
            v = -1i*w*u;                                                    % velocity
            Pxz(iz,ix,:) = -v'*sigT/2;                                    % Poynting vector

            uxz(iz,ix,:) = u;
            sigxz(iz,ix,:) = sig;                                         % stress

        end
    end
    pxz = mean(sigxz(:,:,1:3),3);
end

% return real parts only!
uxz = real(uxz);
pxz = real(pxz);
sigxz = real(sigxz);
Pxz = real(Pxz);
end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
