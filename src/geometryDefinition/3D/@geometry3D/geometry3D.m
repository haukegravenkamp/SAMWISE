classdef geometry3D < geometry
    % general 3D geometries
    % subclasses specify either pre-defined geometries or user-defined mesh

    properties
        eleOrder = 'auto'                                                   % interpolation order of element, chosen automatically by default
        materialNumber = 1                                                  % material number of entire geometry if not specified otherwise
        curvedEdgeFunc                                                      % function defining a curved edge
        scaleGeometry = [1, 1];                                             % scale entire cross-section by these factors in y and z direction
    end

    properties (Hidden = true)
        spatialDimension = 3                                                % number of dimension, relevant for PDE
    end
    
    methods
        function obj = geometry3D
            
        end
        
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
