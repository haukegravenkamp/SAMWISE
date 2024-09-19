clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'steel';                                                % name for reading from database

%% geometry
geom = Lshape;                                                              % create geometry
geom.Ly = 50;
geom.Lz = 50;
geom.hy = 10;
geom.hz = 10;
geom.materialNumber = 1;

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.numerics.eleOrderFunction = @(w) ceil(3 + 0.5*w);

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 0.1;
sol.nSteps = 100;
% sol.nSolutions = 40;

%% call program
[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(materials, geom, opt, sol);

%% plot
figure
plot(sol,'all',opt);

%% compare with old results
ref = load('Lshape_steel.mat');
figure
hold all
plot(ref.f,ref.cp,'.k','LineWidth',1.5)

[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o'});
a = gca;
legend([a.Children(1),a.Children(end)],{'2024','2016'},'Location','best');


