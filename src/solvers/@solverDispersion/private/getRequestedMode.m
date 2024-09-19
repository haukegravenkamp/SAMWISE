

function [uPlot,omegaKmode,indOmK]=getRequestedMode(omega,k,u,omegaPlot,kPlot)
[~,indOmega] = min(abs(omega-omegaPlot));                                   % find frequency closest to requested value
omegaMode = omega(indOmega);                                                % actual frequency
if isreal(kPlot)
    [~,indk] = mink(abs(real(k(:,indOmega))-kPlot),2);                % find two wavenumbers closest to requested value
else
    [~,indk] = mink(abs(k(:,indOmega)-kPlot),2);                            % find two wavenumbers closest to requested value
end
k1 = abs(k(indk(1),indOmega));
k2 = abs(k(indk(2),indOmega));

if abs((k1-k2)/k1) < 1e-6                                                   % if absolute values are practically the same (likely different sign of attenuation)
    [~,indA] = max(imag(k(indk,indOmega)));
    indk = indk(indA);
else
    indk = indk(1);
end
kMode = k(indk,indOmega);                                                   % actual wavenumber
uPlot = u(:,indk,indOmega);                                                 % mode shape
omegaKmode = [omegaMode, kMode];                                            % store omega,k
indOmK = [indOmega, indk];                                                  % store indices
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
