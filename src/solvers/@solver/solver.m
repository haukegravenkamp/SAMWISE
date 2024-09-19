classdef solver < matlab.mixin.Heterogeneous
    % SOLVER generic solver 
    %
    % fMax: maximum frequency, defined in parent class as it is used to
    % determine the element order automatically; hence it should exist in
    % all solvers - at least if automatic element order computation is
    % active.
    properties
        fMax = 'auto'                                                       % maximum temporal frequency (auto: choose such that dimensionless frequency is 10)
        tCPU                                                                % OUTPUT: computing time
    end
    
    methods
        function obj = solver
            
        end
    end
    methods (Static)
        Z = Zmatrix(omega,E0,E11,E12,E2,M0,C0)
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
