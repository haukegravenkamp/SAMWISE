
%% Example 5.1 Homogeneous plate
% Gravenkamp, Hauke, Chongmin Song, and Jens Prager. “A Numerical Approach 
% for the Computation of Dispersion Relations for Plate Structures Using 
% the Scaled Boundary Finite Element Method.” Journal of Sound and 
% Vibration 331 (2012): 2543–57. https://doi.org/10.1016/j.jsv.2012.01.029.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials
materials = material;                                                       % create array of materials
materials.elastic.G   = 1;                                                  % shear modulus
materials.elastic.nu  = 0.3;                                                % Poisson's ratio
materials.elastic.rho = 1;                                                  % density
% See example_plate_Brass.m. Only, here we define the material parameters
% explicitly, as this material is not in the database. Note that you can
% choose almost any valid combination of elastic constants, the remaining
% ones will be set automatically. Instead of (G,nu), you may set 
% (E,nu), (E,G), (G,lambda), (cs,cl), (cs,nu), (cl,nu)
% or even directly set the elasticity matrix in Voigt notation.

%% geometry
geom = plate;                                                               % initialize plate geometry
geom.layers.thickness = 1;                                                  % thickness
geom.layers.material  = 1;                                                  % material number

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 3.2;
sol.nSteps = 50;

%% options
opt = option;
opt.plotting.angularFrequency = true;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, sol);

%% compare with analytical solution
ref = load('plateHomogeneous_nu03.mat');
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1)
axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
legend([axH{1}.Children(1),axH{1}.Children(end)],{'samwise','analytical'},'Location','best');
ylim([0 4.5])





