classdef geometry2D < geometry
    % general 2D geometry, used for plates and cylinders (!)
    
    properties
        layers = layer                                                      % array of layers
        r0 = 0                                                              % vertical position of first layer
    end

    properties (Hidden = true)
        thicknessTotal
    end
    
    methods
        function obj = geometry2D
            
        end
        
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
