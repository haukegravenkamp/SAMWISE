classdef bcUnbounded < bc
    % coupling to unbounded medium

    properties
        material                                                            % material number of unbounded medium
    end
    properties (Hidden)
        dofsUnbounded                                                       % dofs corresponding to unbounded domain (only if additional unknowns are introduced by solver)
        rMatrices                                                           % matrices R0,...,R4 relating displacements and tractions
    end

    methods
        function obj = bcUnbounded
        end
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
