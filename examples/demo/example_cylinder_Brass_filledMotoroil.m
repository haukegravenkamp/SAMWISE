clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'Brass';                                                % name for reading from database
materials(2).name = 'motorOil';                                                % name for reading from database

%% geometry
geom = cylinder;                                         % create geometry of type cylinder

layerThickness = 1;

geom.layers(1).thickness = layerThickness;                                  % thickness
geom.layers(1).material = 2;                                                % material number
geom.layers(2).thickness = layerThickness;                                  % thickness
geom.layers(2).material = 1;                                                % material number
geom.r0 = 0;                                                                % inner radius

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.model.circumferentialOrders = 0:1;

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 2/layerThickness;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
plot(sol,'all',opt);

%% compare with disperse
if layerThickness == 2
    ref = load('cylinder_MotoroilBrass_2mm2mm_disperse_cp.mat');
else
    ref = load('cylinder_MotoroilBrass_ri_0mm_disperse_cp.mat');
end

figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1.5)

hold all
[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
a = gca;
legend([a.Children(1),a.Children(end)],{'samwise','disperse'},'Location','best');


