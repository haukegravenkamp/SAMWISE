classdef pdeAcoustics < pde
    %% pdeAcoustics
    % Partial differential equations defining linear acoustics.
    %
    %   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
    properties (Constant)
        name = 'linear acoustics';                                          % informative name
    end 
    
    properties
        dof = 1;                                                            % one DOF per node
        nTerms = 2;                                                         % number of terms
        bMatrices                                                           % cell with "b-matrices", one term per row
        ordersSpaceTime = [1 1 0; 0 0 2];                                   % dx test functions, dx trial functions, dtt 
        variableName  = '$p$'                                               % name of unknowns, mainly for plotting
        colormap                                                            % default colormap in wavefield plots
    end
    
    methods
        function obj = pdeAcoustics(varargin)                               % constructor
            obj = obj@pde;                                                  % call superclass constructor
        end

        function b = get.bMatrices(obj)
            spDim = obj.spatialDimension;
            % isCylinder = obj.isCylinder;
            switch spDim
                case 2
                    bx = [1 0]';
                    by = [0 0]';
                    bz = [0 1]';
                case 3
                    bx = [1 0 0]';
                    by = [0 1 0]';
                    bz = [0 0 1]';

            end
            b = {bx,by,bz,[],bx,by,bz,[]; [],[],[],[],[],[],[],[]};

        end

        function c = get.colormap(~)
                c = myBlues(256);
        end

    end
    methods(Static)
        
    end
    methods(Sealed)
        
    end
    
end


