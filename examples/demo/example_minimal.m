%% minimal examples
% This is a minimal example in the sense that it only creates the objects
% that are absolutely essential to run the code, everything else is set to
% default values

%% geometry

geom = plate; 
% create geometry of type plate
% 'plate' is a subclass of geometry2D. It defines the basic properties of
% the geometry. To see which other geometries there are, you can check out
% the folder \src\geometryDefinition

geom.layers.thickness = 1;           
% create layer of thickness 1 
% plates and cylinders are defined by layers. The plate is automatically
% created with one layer, more can be added, see the following examples.
% We normally simply use one layer for each material. There is no need to 
% sub-divide a physical layer into several finite elements! 

geom.layers.material = 1;                                                   
% assign material number to layer
% As we will see in more complex examples, we can create a vector of
% different materials. Here, we assign a material to a layer by the number
% of the material in the vector. In this minimal example, we do not
% explicitly define a material, it will default to one with a shear modulus
% of 1, density of 1, and Poisson's ratio of 1/3.

%% solver

sol = solverDispersion;                                                     
% create solver
% Create an instance of the class 'solverDispersion' with default settings.
% By default, solutions are computed at 200 frequencies, between 0 and a
% dimensionless frequency of 10. The dimensionless frequency is defined as 
%      a = 2*pi*f*L/c 
% with f: temporal frequency, L: thickness, c: smallest wave velocity

%% call samwise

[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(geom,sol);               
% call program
% we pass the objects created above to the main program. The arguments can
% be passed in any order. All optional objects (boundary conditions,
% options, materials) are set to default values.

%% plot
plot(sol);                                                                  
% plot solution
% In our case, this will call the function 'plotSolver' of the
% solverDispersion class. It plots the dispersion curves in terms of
% wavenumber, phase velocity, group velocity, attenuation





