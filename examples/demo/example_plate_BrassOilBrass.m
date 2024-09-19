clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'Brass';                                                % name for reading from database
materials(2).name = 'motorOil';                                                % name for reading from database

%% geometry
geom = plate;                                            % create geometry of type cylinder

geom.layers(1).thickness = 1;                                               % thickness
geom.layers(1).material = 1;                                                % material number

geom.layers(2).thickness = 1;                                               % thickness
geom.layers(2).material = 2;                                                % material number

geom.layers(3).thickness = 1;                                               % thickness
geom.layers(3).material = 1;                                                % material number

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.numerics.decomposeEVP = false;

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 2.4;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
plot(sol,'all',opt);

%% compare with disperse
ref = load('layers_BrassOilBrass.mat');
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1.5)
hold all
[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
a = gca;
legend([a.Children(1),a.Children(end)],{'samwise','disperse'},'Location','best');


