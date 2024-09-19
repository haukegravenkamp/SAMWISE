classdef optionNumerics < handle & matlab.mixin.SetGet
    % optionsNumerics
    % settings regarding general numerical properties

    properties
        shpFunType = 'Lagrange'                                             % type of shape function
        intType = 'auto'                                                    % type of integration; 0: Gauss, 1: Gauss-Lobatto
        addIntPoints = 0                                                    % additional integration points per direction
        lumpingM = false                                                    % lump mass matrix
        lumpingE0 = false                                                   % lump E0 matrix
        eleOrderFunction = @(w) ceil(3 + 0.5*w)                             % function to set element order depending on dimensionless frequency w
        eleOrderDefault = 5                                                 % if everything else fails, use this element order
        hRefinement = 0                                                     % number of refinement steps, each step reduces the element length by 1/2
        computeEigenvectors = true                                          % whether to compute eigenvectors (mode shapes), may increase computational costs
        decomposeEVP = true                                                 % try to decompose quadratic eigenvalue problem for more efficient linearization
        unboundedRemoveIncoming = true                                      % when coupling with unbounded domains, remove modes propagating towards the waveguide
        unboundedUsePoyntingVec = true;                                     % when removing incoming waves, use  Poynting vector rather than just sign of wavenumber
        unboundedRemoveThreshold = 1e-3;                                    % threshold for determining direction of modes, relative to total maximum                                            
        removeNegativeAttenuation = true;                                   % whether to remove modes with negative attenuation after solution
        computeResidual = false;                                            % whether to compute residuals of solution (if supported by solver)
        nEigenvaluesDefault = 50;                                           % if not specified, compute this many eigenvalues (in case of large, sparse matrices)
        useSparseSize = 200;                                                % use sparse matrices and 'eigs' for matrices above this size
    end

    methods
        function obj = optionNumerics

        end

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
