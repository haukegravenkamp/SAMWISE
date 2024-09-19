classdef pdeElasticity < pde
    %% pdeElasticity
    % Partial differential equations defining linear elasticity.
    %
    %   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com

    properties (Constant)
        name = 'Linear Elasticity';                                         % name of pde
    end

    properties
        dof = 3;                                                            % specified in 'definePde'
        nTerms = 2;                                                         % number of terms
        bMatrices                                                           % cell with "b-matrices", one term per row,
        assumption2D = 'none'                                               % behavior of 2D bodies
        ordersSpaceTime = [1 1 0; 0 0 2];                                   % dx test functions, dx trial functions, dtt
        variableName = '$u$'                                                % name of unknowns, mainly for plotting
        colormap                                                            % default colormap in wavefield plots
        plotDof = [1 3];                                                    % default dofs for plotting wave fields
    end

    methods
        function obj = pdeElasticity(assumption2D)                          % constructor
            obj=obj@pde;                                                    % call superclass constructor
            if nargin > 0
                obj.assumption2D = assumption2D;
                if strcmp(assumption2D,'planeStrain') || ...
                        strcmp(assumption2D,'planeStress')
                    obj.dof = 2;
                elseif strcmp(assumption2D,'out-of-plane')
                    obj.dof = 1;
                    obj.plotDof = 1;
                end
            end
        end


        function b = get.bMatrices(obj)
            spDim = obj.spatialDimension;                                   % number of spatial dimensions given by geometry
            isCylinder = obj.isCylinder;

            switch spDim
                case 2

                    switch obj.dof                                          % number of dofs, depending on strain state in 2D
                        case 1                                              % out of plane, same as acoustics
                            bx = [1 0]';
                            by = [0 0]';
                            bz = [0 1]';
                        case 2                                              % plane strain or plane stress
                            bx=[1 0; 0 0; 0 1];
                            bz=[0 0; 0 1; 1 0];
                            by = bz *0;
                        case 3                                              % all dofs
                            bx=[1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1; 0 1 0];
                            bz=[0 0 0; 0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 0];
                            by = bz *0;
                    end

                case 3                                                      % 3D models
                    bx = [1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1; 0 1 0];
                    by = [0 0 0; 0 1 0; 0 0 0; 0 0 1; 0 0 0; 1 0 0];
                    bz = [0 0 0; 0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 0];

            end
            b = {bx,by,bz,[],bx,by,bz,[]; [],[],[],[],[],[],[],[]};

            if isCylinder
                br = [0 0 0; 0 0 1; 0 0 0; 0 -1 0; 0 0 0; 0 0 0];
                b{1,4} = br;
                b{1,8} = br;
            end


        end

        function c = get.colormap(~)
            % c = myGolds(256);
            % c = myGreens(256);
            % c = myBlues(256,[2,90]);
            c=gray;
        end



    end
    methods(Static)

    end
    methods(Sealed)

    end

end

 
