classdef cylinder < geometry2D
    % cylinder geometry, consists of layers
    % inherits from plates, as the definition of layers is the same

    properties
        assumption2D = 'cylinder'                                        % 'none', 'planeStrain', 'planeStress' (if applicable)
        spatialDimension = 3                                                % number of dimension, relevant for PDE
        % spatialDimension is 2 for plate, even though out-of-plane motion
        % may be included, as there is no derivative considered

    end

    properties (Hidden = true)
        radiusOuter                                                         % outer radius (radiusInner + thicknessTotal)
    end

    methods
        function obj = cylinder
            
        end

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
