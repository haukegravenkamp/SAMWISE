
%% geometry
geom = plate;                                            % create geometry of type plate
geom.layers.thickness = 1;                                                  % create layer of thickness 1 
geom.layers.material = 1;                                                   % assign material number to layer

%% solver
sol = solverDispersion;                                                     % create solver
sol.fMin = 0;                                                               % minimum frequency

[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(geom,sol);               % call program

plot(sol);                                                                  % plot solution





