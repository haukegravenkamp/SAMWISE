clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % name for reading from database

%% geometry
geom = cylinder;

geom.layers.thickness = 1;                                                  % thickness
geom.layers.material = 1;                                                   % material number
geom.r0 = 0;                                                                % inner radius, 0 or 1

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.model.circumferentialOrders = 0:1;

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 7;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
plot(sol,'all',opt);

%% compare with disperse
ref = load(['cylinder_Brass_ri_',num2str(geom.r0),'mm_disperse_cp.mat']);

figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1.5)

[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
a = gca;
legend([a.Children(1),a.Children(end)],{'samwise','disperse'},'Location','best');


