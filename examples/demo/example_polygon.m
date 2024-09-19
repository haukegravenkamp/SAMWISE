clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % name for reading from database

%% geometry
geom = polygonRegular;                                                      % regular polyon
geom.R = 1;                                                                 % outer radius
geom.n = 7;                                                                 % number of corners

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.numerics.eleOrderFunction = @(w)ceil(3+0.5*w);
opt.numerics.nEigenvaluesDefault = 30;

%% solver
sol = solverDispersion;
sol.fMin = 1e-3;
sol.fMax = 0.8;
sol.nSteps = 50;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
figure
plot(sol,'all',opt);
drawnow

%% compare with cylinder computed using disperse
% agreement is better the more corners the polygon has

ref = load('cylinder_Brass_ri_0mm_disperse_cp.mat');
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1.5)
[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
a = gca;
legend([a.Children(1),a.Children(end)],{'samwise','disperse'},'Location','best');
xlim([sol.fMin,sol.fMax])

