
%% Example 5.2 Layered composite
% Gravenkamp, Hauke, Chongmin Song, and Jens Prager. “A Numerical Approach 
% for the Computation of Dispersion Relations for Plate Structures Using 
% the Scaled Boundary Finite Element Method.” Journal of Sound and 
% Vibration 331 (2012): 2543–57. https://doi.org/10.1016/j.jsv.2012.01.029.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % material from database
materials(2).name = 'titanium';                                             % material from database

%% geometry
geom = plate;                                                               % initialize plate geometry
geom.layers(1).thickness = 1/5;                                             % thickness
geom.layers(1).material  = 1;                                               % brass
geom.layers(2).thickness = 1/5;                                             % thickness
geom.layers(2).material  = 2;                                               % titanium
geom.layers(3).thickness = 1/5;                                             % thickness
geom.layers(3).material  = 1;                                               % brass
geom.layers(4).thickness = 1/5;                                             % thickness
geom.layers(4).material  = 2;                                               % titanium
geom.layers(5).thickness = 1/5;                                             % thickness
geom.layers(5).material  = 1;                                               % brass
% Here, the plate consists of five layers; hence, we keep extending the
% layer object and assign the materials 1 and 2 in turn. Each layer has a
% thickness of 1/5 mm, so that the total thickness of the plate is 1 mm.


%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 10.5;
sol.nSteps = 50;

%% options
opt = option;
opt.plotting.angularFrequency = ~true;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, sol);

%% compare with analytical solution
ref = load('composite_BrTiBrTiBr.mat');
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1)
axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
legend([axH{1}.Children(1),axH{1}.Children(end)],{'samwise','disperse'},'Location','best');
xlim([sol.fMin,sol.fMax])
ylim([0 11])





