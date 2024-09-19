classdef solverDispersion_embeddedExact < solverDispersion
    % solverDispersion_embeddedExact
    % solver for computing dispersion curves of plate structures coupled to
    % infinite media
    % 

    properties 
        kz_top                                                              % vertical wavenumber in unbounded domain, coupled to top surface
        kz_bottom                                                           % vertical wavenumber in unbounded domain, coupled to bottom surface
        Pz_top                                                              % power flux in vertical direction, top surface
        Pz_bottom                                                           % power flux in vertical direction, bottom surface
    end

    properties %(Access = protected) 
        Rp                                                                  % coupling matrices, corresponding to pressure waves
        Rs                                                                  % coupling matrices, corresponding to shear waves
        waveVel                                                             % wave velocities in unbounded medium [bottom:cl, bottom:cs; top:cl, top:cs]
        sameUnbMedia = false                                                % whether plate is coupled to the same unbounded medium at both surfaces
    end
    


    methods
        function obj = solverDispersion_embeddedExact

        end

    end
    methods (Static)
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
