classdef material 
    % MATERIAL 
    
    properties
        name                                                                % arbitrary name (string)
        behavior = 'elastic'                                                % 'elastic', 'acoustic' (i.e. one of the properties listed below)
        elastic = parametersElastic                                         % object containing elastic constants
        acoustic = parametersAcoustic                                       % object containing acoustc constants
        
    end

    properties (Hidden = true)
        parameters                                                          % relevant parameters depending on material behavior
    end
    
    methods
        function obj = material

        end
    end
    methods (Static)
        databasePath = getDatabasePath;
    end

end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
