classdef optionPlotting < handle & matlab.mixin.SetGet

    properties
        dimless = false                                                     % whether to use dimensionless properties for plotting
        angularFrequency = false                                            % whether to use dimensionless frequency for plotting
        maxAttenuation = 0;                                                 % plot only modes with attenuation smaller than this value; 0 for propagating modes (default)                                                      
        propagatingModesThreshold = 1e-6;                                   % attenuation threshold for modes being considered as propagating
        omegaVSk = false                                                    % whether to plot frequency vs wavenumber (true) or vice versa (false)
        plotMesh = 1                                                        % whether to plot the mesh
        plotNodes = 'auto';                                                 % in mesh plots, whether to plot nodes (auto: only if less than 500)
        fillElements = 0                                                    % whether to fill the elements in mesh with color
        wavefieldPointsPerWavelength = 20                                   % number of points per smallest wavelength in wavefield plots
        waveFieldLx = 2                                                     % size of plotting domain in x-direction, relative to waveguide's thickness
        waveFieldLz = 1                                                     % size of plotting domain in z-direction in unbounded domains, relative to waveguide's thickness
        waveFieldLineStyle = '-'                                            % line style for surf plot
        suppressOutput = 0                                                  % whether to print data and cpu times on screen
        gifColorbar = false                                                 % include the colorbar when creating gifs using saveGif
        gifMaximize = false                                                 % maximize figure to whole screen when creating gifs (higher resolution, larger file)
        gifBorder = 20                                                      % keep so many white pixels around content when creating gifs
        plotFontSizes = [12,14]                                             % font sizes for ticks and labels
        distortedGrid = true                                                % in wavefield plots, show distortion of plot due to displacements
        distortionFactor = 1                                                % normalize distortion, relative to grid size
        acousticDisplacement = false                                        % in acoustic problem, plot displacement instead of pressure if possible
        plotZeroGraphs = false;                                             % whether to include graphs that only contain zeros
        omitNegativeWavenumber = ~true;                                     % whether to remove modes with negative wavenumber (including backward waves)
        omitNegativeAttenuation = true;                                     % whether to remove modes with negative attenuation

        labelFrequency = '$f$'                                              % standard axis labels to be used in plots 
        labelCircularFrequency = '$\omega$'
        labelWavenumber = '$k$'
        labelPhaseVelocity = '$c_p$'
        labelGroupVelocity = '$c_g$'
        labelAttenuation = '$\eta$'
        labelDisplacement = '$u$'
    end

    methods
        function obj = optionPlotting
        end

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
