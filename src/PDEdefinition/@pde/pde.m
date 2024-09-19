classdef pde < handle & matlab.mixin.Heterogeneous & matlab.mixin.SetGet
    %% PDE
    % class defining partial differential equations.
    % The different pde's are defined as subclasses

    properties (Abstract)
        % NOTE: The ABSTRACT properties must be defined and set by each
        % subclass. In the subclass, they are defined as normal properties
        dof                                                                 % number of dofs per node
        nTerms                                                              % number of terms
        bMatrices                                                           % cell with "b-matrices", one term per row,
        % each row contains {bx1,by1,bz1, bx2,by2,bz2}, where 1,2 are for
        % test and trial functions, respectively
        ordersSpaceTime                                                     % order of terms in frequency
        variableName                                                        % name of unknowns, mainly for plotting

    end

    properties (Abstract, Constant)
        name                                                                % informative name of the problem
    end

    properties
        spatialDimension                                                    % just copied from geometry
        isCylinder = false                                                  % if geometry is cylinder
    end

    properties (SetAccess = protected, Hidden=true)

    end
    properties (Hidden=true)

    end
    methods
        function obj=pde                                                    % constructor

        end
 
    end

    methods(Abstract)

    end
    methods(Static)
        % assignPDE(msh,mat)
        pdes = getAllPdes(geo)

    end
    methods(Sealed)

    end

end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
