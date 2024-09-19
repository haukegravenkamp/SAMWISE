classdef circularHollow < geometry3D
    % pre-defined geometry of a hollow cylinder (pipe) of outer radius Ro
    % and inner Radius Ri
    % This geometry creates a fully discretized cross-section, which is
    % only required if either the material or the boundary condition are
    % not of rotational symmetry. 
    % In all other cases, use the geometry "cylinder", which discretizes
    % only the thickness.

    properties
        Ri
        Ro
    end

    properties (Hidden = true)
        
    end
    
    methods

        
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
