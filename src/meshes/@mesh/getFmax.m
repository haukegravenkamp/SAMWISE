function fMax = getFmax(fMax,eleSizes,cMats,materialNumber)

if isscalar(fMax) && fMax>=0                                                % maximum frequency provided properly
    return
else
    if  ~ischar(fMax)||~strcmp(fMax,'auto')                                 % neither scalar nor 'auto' -> unknown
        warning('solver property fMax not defined properly. Using default (dimensionless frequency aMax = 10)')
    end
    aMax = 10;                                                              % desired dimensionless frequency
    L = eleSizes;                                                           % element size
    cMin = min(cMats,[],2);                                                 % (smallest) velocity of each material
    c = cMin(materialNumber);                                               % (smallest) velocity of each element
    fMax = min(1./(2*pi*L./c/aMax));                                        % adequate frequency for each element
end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
