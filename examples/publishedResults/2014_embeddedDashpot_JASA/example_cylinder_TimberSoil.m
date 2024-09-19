
clc
clear 
close all

%% Example III.C Timber pole embedded in soil
% Gravenkamp, Hauke, Carolin Birk, and Chongmin Song. “Computation of 
% Dispersion Curves for Embedded Waveguides Using a Dashpot Boundary 
% Condition.” The Journal of the Acoustical Society of America 135, 
% no. 3 (2014): 1127–38. https://doi.org/10.1121/1.4864303.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials

mat = material;                                                             % create array of materials
mat(1).name = 'timber';                                                     % name 
mat(1).elastic.isotropy = 2;                                                % orthotropic
mat(1).elastic.E1   = 23;
mat(1).elastic.E2   = 1.177;
mat(1).elastic.E3   = 2.665;
mat(1).elastic.G12  = 1.403;
mat(1).elastic.G23  = 0.483;
mat(1).elastic.G31  = 2.047;
mat(1).elastic.nu12 = 0.43;
mat(1).elastic.nu13 = 0.35;
mat(1).elastic.nu23 = 0.309;
mat(1).elastic.rho  = 0.9;

mat(2).name = 'soil';           
mat(2).elastic.G   = 0.04;
mat(2).elastic.rho = 1.5;
mat(2).elastic.nu  = 0.3;

%% geometry

geom = cylinder;                                                            % create geometry of type cylinder
geom.layers.thickness = 150;                                                % thickness
geom.layers.material  = 1;                                                  % material number of the layer

%% boundary conditions

bCond = bcUnbounded;                                                        % coupling to unbounded domain
bCond.location = 'top';                                                     % upper surface
bCond.material = 2;                                                         % material number of unbounded domain (refers to materials array above)

%% settings 

opt = option;
opt.model.unboundedModel = 'dashpot';                                       % approximate dashpot boundary conditions -> standard solver
opt.plotting.maxAttenuation = 120;                                          % maximum attenuation, just for filtering/plotting

%% solver

sol = solverDispersion;
sol.fMin = 0.0;
sol.fMax = 0.01;
sol.nSteps = 200;

%% call program

[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(mat, geom, opt, bCond, sol);

%% compare with absorbing region
ref = load('cylinder_timberSoil.mat');
figure
hold all
axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
plot(ref.fcp/1e3,ref.cp,'xk','LineWidth',1)
legend([axH{1}.Children(end),axH{1}.Children(1)],{'samwise','absorbing region'},'Location','best');
ylim([0 6])
xlim([0 0.01])

figure
hold all
axH = plot(sol,'att',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
plot(ref.fatt/1e3,ref.att,'xk','LineWidth',1)
legend([axH{1}.Children(end),axH{1}.Children(1)],{'samwise','absorbing region'},'Location','best');
xlim([0 0.01])


