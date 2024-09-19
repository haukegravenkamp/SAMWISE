function [eta,wgt]=readIntegrationPoints(obj,nP,intType)

dimension=numel(nP);                                                        % number of points MUST be provided for each direction in a vector.
                                                                            % If it's a scalar, we assume a 1D problem.
switch dimension
    case 1                                                                  % 1D
        [eta,wgt]=intPoints1D(obj,nP,intType);
    case 2                                                                  % 2D
        [eta,wgt]=intPoints2D(obj,nP(1),nP(2),intType);
    otherwise
        error('function readIntegrationPoints not implemented for this case')
end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
