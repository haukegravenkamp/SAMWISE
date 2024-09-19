classdef optionModel < handle & matlab.mixin.SetGet
    % optionsModel
    % settings regarding the physical model
    % 

    properties
        plateModeTypes = 'auto'                                             % FOR PLATES, which modes to compute: 0 all modes; 1 Lamb modes; 2 SH (if applicable)
        cylinderModeTypes = 'auto'                                          % FOR CYLINDERS, which modes to compute: 0 all modes; 1 L,F; 2 torsional (if applicable)
        circumferentialOrders = 0:1                                      % FOR CYLINDERS: order(s) in circumferential direction, 1D array
        unboundedModel = 'dashpot'                                          % how to describe unbounded domain coupled to waveguide's surface: 'dashpot', 'exact'
        rayleighDamping = {[],[],[],[],[],[],[],false}                      % parameters for Rayleigh damping
    end

    methods
        function obj = optionModel

        end

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
