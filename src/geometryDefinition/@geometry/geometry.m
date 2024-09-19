classdef geometry
    % geometry
    % general geometry class
    % constructor calls subclasses based on geometry type (plate, cylinder, 3D)

    properties
        
    end

    properties (Abstract)
        spatialDimension

    end


    methods
        % constructor
        function obj = geometry

        end

    end

    methods(Static)
        % creating different types of geometry
        function obj = create(geometryType)
            switch geometryType
                case 'plate'
                    obj=plate;
                case 'cylinder'
                    obj=cylinder;
                case 'none'

            end

        end
    end



end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
