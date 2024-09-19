function [uxz,pxz,sigxz,Pxz] = pointValuesUnbounded(obj,uEle,omegaKmode,indOmK,x,z,nx,nz,matEle,unbBTa,pdeName)

indOmega = indOmK(1);
indk = indOmK(2);

w = omegaKmode(1);
kx = omegaKmode(2);

dof = numel(uEle);
uxz = zeros(nz,nx,3);                                                       % allocate displacement
Pxz = zeros(nz,nx,3);                                                       % allocate Poynting vector

rho = matEle.parameters.rho;                                                % read mass density

nU = numel(uEle);

%% check properties of unbounded domain
if unbBTa == 1                                                              % bottom
    kzAll = obj.kz_bottom;                                                  % read wavenumbers in unbounded domain
    z0 = max(z);                                                            % interface position
elseif unbBTa == 2                                                          % top
    kzAll = obj.kz_top;
    z0 = min(z);
end
if ismatrix(kzAll)                                                          % only one wavenumber per mode and frequency
    kz = squeeze(kzAll(indk,indOmega,:));
elseif ndims(kzAll) == 3                                                    % two wavenumbers
    kz = squeeze(kzAll(indk,indOmega,:));
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if strcmp(pdeName,'pdeAcoustics')
else
end


if strcmp(pdeName,'pdeAcoustics')
    % only for acoustic problem: particle acceleration is related
    % to pressure:
    % a = -1/rho grad(p) -> u = 1/(omega²*rho) grad(p)
    pxz = zeros(nz,nx);                                                     % allocate acoustic pressure
    for iz = 1:nz                                                           % loop vertical position in element
        for ix = 1:nx                                                       % loop horizontal position
            p = uEle * exp(1i*kx*x(ix)) * exp(1i*kz*(z(iz)-z0));            % solution (pressure)
            pdx =  p*1i*kx;                                                 % dp/dx
            pdz =  p*1i*kz;                                                 % dp/dz
            gradP = [pdx,0,pdz];                                            % grad(p)
            u = 1/w^2/rho * gradP;                                          % displacement
            v = -1i*w*u;                                                    % velocity
            
            Pxz(iz,ix,:) = real(v'*p)/2;                                   % Poynting vector
            pxz(iz,ix)   = p;
            uxz(iz,ix,:) = u;
        end
    end
    sigxz = [];
end

if strcmp(pdeName,'pdeElasticity')
    D = matEle.elastic.D;                                                   % read elasticity matrix
    Dk = diag(kz([1,2,2]));

    switch dof
        case 1                                                              % out of plane
            uEle = reshape([zeros(2,nU);uEle.'],[],1);                      % only third component is nonzero
        case 2                                                              % plane strain/stress
            unonzero = reshape(uEle,2,[]);
            uEle = reshape([unonzero(1,:); zeros(1,nU/2); unonzero(2,:)],[],1);       % components 1 and 3 nonzero
            % uEle = reshape([reshape(uEle,2,[]); zeros(1,nU/2)],[],1);       % components 1 and 2 nonzero
    end

    for iz = 1:nz                                                           % loop vertical position in element
        for ix = 1:nx                                                       % loop horizontal position
            A = 1i*kx*diag([1 1 -1]) + 1i*[0 0 kz(2); 0 0 0; kz(1) 0 0];
            B = diag(exp(1i*(z(iz)-z0)*kz([1,2,2])));
            u = A * B * exp(1i*kx*x(ix)) * uEle;
            udx = 1i*kx*u;
            udz = A * 1i*Dk * B * exp(1i*kx*x(ix)) * uEle;
            udy = 0*u;
            e = [udx(1), udy(2), udz(3),...
                udy(3)+udz(2), udx(3)+udz(1), udy(1)+udx(2)].';             % strain
            sig = D*e;                                                      % stress (Voigt notation)
            sigT = voigt(sig);                                              % stress tensor

            v = -1i*w*u;                                                    % velocity
            Pxz(iz,ix,:) = -v'*sigT/2;                                      % Poynting vector

            uxz(iz,ix,:) = u;
            sigxz(iz,ix,:) = sig;                                           % stress

        end
    end
    pxz = mean(sigxz(:,:,1:3),3);
end

% return real parts
uxz = real(uxz);
pxz = real(pxz);
sigxz = real(sigxz);
Pxz = real(Pxz);


end
