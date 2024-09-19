classdef plate < geometry2D
    % plate geometry, consists of layers

    
    properties
        assumption2D = 'planeStrain'                                        % 'none', 'planeStrain', 'planeStress' (if applicable)
        spatialDimension = 2                                                % number of dimension, relevant for PDE
        % spatialDimension is 2 for plate, even though out-of-plane motion
        % may be included, as there is no derivative considered
    end

    properties (Hidden = true)
        
    end
    
    methods
        function obj = plate
            
        end
        
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
