classdef glbSys
    % glbSys
    % storing the element matrices and related properties 
    % 
    
    properties
        feMatrices ={}                                                      % cell containing all finite-element matrices                                              
        fixedDofs                                                           % indices of dofs with vanishing solution
        bcDof                                                               % indices associated with each boundary condition (logical)
        shpDy                                                               % derivative of shape function in y-direction at nodes
    end
    
    methods
        function obj = glbSys()

        end
        
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
