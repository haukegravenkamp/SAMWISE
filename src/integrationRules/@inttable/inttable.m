classdef inttable
    %% class for integration table
    % contains the weights and integration points for Gauss and GLL
    % quadrature.
    %
    %   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
    
    properties
        order
        Gauss
        GLL
    end
    
    methods
        function obj=inttable(order) % constructor
            
            if nargin>0
                nPoints = order+1;
                % Gauss Quadrature
                for i = nPoints:-1:1
                    obj(i).order = i-1;
                    obj(i).Gauss = inttable.Gpoints(i);
                end
                
                % Gauss-Lobatto Quadrature
                for i = nPoints:-1:1
                    obj(i).order = i-1;
                    obj(i).GLL = inttable.GLLpoints(i);
                    
                end
                
            end
        end
    end
    
    methods(Static)
        Gauss=Gpoints(N);
        GaussLobatto=GLLpoints(N) ;
    end
end
 
