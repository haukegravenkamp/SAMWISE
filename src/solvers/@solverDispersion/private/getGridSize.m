
function gridSize = getGridSize(obj,kMode,indOmK,wavefieldPointsPerWavelength,unbBTa)

if (nargin < 5) || isempty(unbBTa)
    unbBTa = false;
end

indOmega = indOmK(1);
indk = indOmK(2);
kMax = abs(real(kMode));
if any(unbBTa==1)                                                           % coupling at bottom
    if ndims(obj.kz_bottom)==3
        kzBot = max(abs(real(obj.kz_bottom(indk,indOmega,:))));             % read wavenumbers in unbounded domain
    elseif ismatrix(obj.kz_bottom)
        kzBot = max(abs(real(obj.kz_bottom(indk,indOmega))));               % read wavenumbers in unbounded domain
    end
    kMax = max(abs(real([kzBot,kMax])));
end
if any(unbBTa==2)                                                           % top
    if ndims(obj.kz_top)==3
        kzTop = max(abs(real(obj.kz_top(indk,indOmega,:))));                % read wavenumbers in unbounded domain
    elseif ismatrix(obj.kz_top)
        kzTop = max(abs(real(obj.kz_top(indk,indOmega))));                  % read wavenumbers in unbounded domain
    end
    kMax = max(abs(real([kzTop,kMax])));
end
gridSize = 2*pi/kMax/wavefieldPointsPerWavelength;


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
