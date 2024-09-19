clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'BrassCui';                                            
cs = 3.2;
materials(1).elastic.cl = 5.9;
materials(1).elastic.cs = cs;
materials(1).elastic.rho = 7.9;

materials(2).name = 'waterCui';                                            
materials(2).behavior = 'acoustic';
rho = 1;
cl = 1.5;
materials(2).acoustic.rho = rho;
materials(2).acoustic.cl = cl;
materials(2).acoustic.K = cl*cl*rho;

%% geometry
geom = cylinder;                                         % create geometry of type cylinder

geom.layers(1).thickness = 9.5;                                  % thickness
geom.layers(1).material = 2;                                                % material number

d = 0.5;
geom.layers(2).thickness = d;                                  % thickness
geom.layers(2).material = 1;                                                % material number
geom.r0 = 0;                                                                % inner radius

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.model.circumferentialOrders = 0;
opt.numerics.eleOrderFunction = @(w)ceil(3+0.6*w);
opt.plotting.omegaVSk = true;
opt.plotting.maxAttenuation = 10;
%% solver
sol = solverDispersion;
sol.fMin = 5.4*cs/d/2/pi;
sol.fMax = 5.6*cs/d/2/pi;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
plot(sol,'all',opt);

%% compare with literature
ref = load('cylinder_waterFilled.mat');
figure
hold all
plot(ref.k,ref.w,'.k')
yl = ylim;
xl = xlim;
% non-dimensionalize
sol.f = 2*pi*sol.f/cs*d;
sol.k = sol.k*d/2/pi;
axH = plot(sol,'k',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
legend([axH{1}.Children(1),axH{1}.Children(end)],{'samwise','Cui et al.'},'Location','best');
xlim([0 0.4])
ylim([5.4 5.6])
xlabel('$k h / 2 \pi$')
ylabel('$\omega h / c_t$')
