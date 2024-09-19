classdef option < handle & matlab.mixin.SetGet
       
    properties

        % model properties
        model = optionModel
        
        % numerical properties
        numerics = optionNumerics

        % post processing
        postprocessing = optionPostprocessing

        % plotting
        plotting = optionPlotting;
        
    end

    properties(Constant)

    end
    
    methods

    end
    
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
