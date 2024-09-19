classdef Lshape < geometry3D
    % pre-defined geometry of an L-shape, characterized by the total
    % dimensions along the two principal axes and the thicknesses hy,hz as
    % shown below:
    % 
    %       Ly
    %     _______
    % hz |_____  |
    %          | | Lz
    %          | |
    %          |_|
    %          hy


    properties
        Ly
        Lz
        hy
        hz
    end

    properties (Hidden = true)
        
    end
    
    methods
        function obj = Lshape
            
        end
        
    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
