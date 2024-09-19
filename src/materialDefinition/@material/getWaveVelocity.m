function cMats = getWaveVelocity(obj,usedMaterials)

nMat = numel(obj);                                                          % number of materials

cMats = zeros(nMat,1);                                                      % relevant wave velocity of each material
for iMat = 1:nMat                                                           % loop materials

    parameters = obj(iMat).parameters;                                      % get material parameters
    if isprop(parameters,'cs') && isfloat(parameters.cs) ...
            && ~isempty(parameters.cs)                                      % shear velocity defined
        cMats(iMat) = parameters.cs;
    elseif isprop(parameters,'cl') && isfloat(parameters.cl) ...
            && ~isempty(parameters.cl)                                      % longitudinal velocity defined
        cMats(iMat) = parameters.cl;
    else                                                                    % neither wave velocities defined
        cMats(iMat) = nan;                                                  % will resort to default later
        if ismember(iMat,usedMaterials)                                     % if material is actually being used
            warning(['material ', num2str(iMat), ' does not define a ' ...   % warn user
                'wave velocity assigned for computing the ' ...
                'dimensionless frequency'])
        end
    end

end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
