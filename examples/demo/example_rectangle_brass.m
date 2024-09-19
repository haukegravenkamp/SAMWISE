clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % name for reading from database

%% geometry
geom = rectangular;                                                              % create geometry
geom.Ly = 1;
geom.Lz = 1;
geom.materialNumber = 1;

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 3;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, sol,opt);

%% plot
plot(sol,'all',opt);

%% compare with 2016 solution
ref = load('square_brass.mat');
figure
hold all
plot(ref.f,ref.cp,'.k','LineWidth',1.5)
hold all
[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o'});
a = gca;
legend([a.Children(1),a.Children(end)],{'2024','2016'},'Location','best');


