classdef layer
    % LAYER Definition of a layer in a plate structure or cylinder


    properties
        thickness = 1
        material = 1                                                        % material number, refers to array of material objects
        eleOrder = 'auto'                                                   % interpolation order of element, chosen automatically by default
        nEle = 1;                                                           % number of elements in layer (typically just one with automatic order)
    end

    methods
        function obj = layer(thickness, material)
            % constructor
            if nargin > 0
                obj.thickness = thickness;
            end
            if nargin > 1
                obj.material = material;
            end
        end

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
