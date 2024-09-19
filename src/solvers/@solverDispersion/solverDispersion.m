classdef solverDispersion < solver
    % solverDispersion
    % parent solver for computing dispersion curves
    % When using this general solver, the adequate child class is chosen
    % automatically. Alternatively, the children can be used directly.

    properties 
        fMin = 0                                                            % INPUT: minimum temporal frequency
        nSteps = 200                                                        % INPUT: number of frequency steps
        nSolutions = 'auto'                                                 % INPUT: number of eigenvalues to compute, 'auto', 'all' or integer

        f                                                                   % OUTPUT: array of temporal frequencies
        omega                                                               % OUTPUT: array of circular frequencies
        k                                                                   % OUTPUT: wave numbers
        cp                                                                  % OUTPUT: phase velocities
        cg                                                                  % OUTPUT: group velocities
        att                                                                 % OUTPUT: attenuations
        phi                                                                 % OUTPUT: eigenvectors
        u                                                                   % OUTPUT: eigenvectors, primary variable
        q                                                                   % OUTPUT: eigenvectors, forces
        uDy                                                                 % OUTPUT: eigenvectors, derivative wrt y
        res                                                                 % OUTPUT: residuals

        
    end

    properties (Access = protected) 
        
    end
    


    methods
        function obj = solverDispersion

        end

    end

    methods (Static)
        [ind,nBl] = blockDecomposition(B,thB)
        vgs = computeGroupVelocity(ki,vu,vq,M0,omega);
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
