classdef result < matlab.mixin.Heterogeneous & handle & matlab.mixin.SetGet
    % result

    properties
        frequency
        angularFrequency
        wavenumber
        phaseVelocity
        groupVelocity
        attenuation
    end

    % synonyms of results variables for convenience
    properties (Hidden = true)
        f
        omega
        k
        cp
        cg
        att
    end

    methods
        function obj = result

        end

    end
    methods 
        function k = get.k(obj)
            k = obj.wavenumber;
        end
        function f = get.f(obj)
            f = obj.frequency;
        end
        function omega = get.omega(obj)
            omega = obj.angularFrequency;
        end
        function cp = get.cp(obj)
            cp = obj.phaseVelocity;
        end
        function cg = get.cg(obj)
            cg = obj.groupVelocity;
        end
        function att = get.att(obj)
            att = obj.attenuation;
        end

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
