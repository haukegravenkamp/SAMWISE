clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % name for reading from database

%% geometry
geom = circular;                                                            % full cylinder
geom.R = 1;                                                                 % radius

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.numerics.eleOrderFunction = @(w)ceil(3+0.5*w);

%% solver
sol = solverDispersion;
sol.fMin = 1e-3;
sol.fMax = 1.5;
sol.nSteps = 50;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
plot(sol,'all',opt);

%% compare with disperse
ref = load('cylinder_Brass_ri_0mm_disperse_cp.mat');
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1.5)

[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
a = gca;
legend([a.Children(1),a.Children(end)],{'samwise','disperse'},'Location','best');
xlim([sol.fMin,sol.fMax])


