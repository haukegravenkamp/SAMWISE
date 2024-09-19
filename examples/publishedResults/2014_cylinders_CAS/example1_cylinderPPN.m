clc
clear 
close all

%% Example 5.1 Natural polypropylene (PPN)
% Gravenkamp, Hauke, Fabian Bause, and Chongmin Song. “On the Computation 
% of Dispersion Curves for Axisymmetric Elastic Waveguides Using the Scaled
% Boundary Finite Element Method.” Computers & Structures 131 (2014): 46–55.
% https://doi.org/10.1016/j.compstruc.2013.10.014.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials
mat = material;                                                             % create array of materials
mat.elastic.cl   = 2.7;                                                     % shear modulus
mat.elastic.nu  = 0.35;                                                     % Poisson's ratio
mat.elastic.rho = 1;                                                        % density

%% geometry
geom = cylinder;                                                            % initialize plate geometry
geom.layers.thickness = 6;                                                  % thickness
geom.layers.material  = 1;                                                  % material number
geom.r0 = 3;

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 1.9;
sol.nSteps = 100;

%% options
opt = option;
opt.model.circumferentialOrders = 0;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(mat, geom, sol,opt);

%% compare with Global Matrix Method
ref = load('cylinderPPN.mat');
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1)
[axH, plH] = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
a = gca;
legend([a.Children(1),a.Children(end)],{'samwise','GMM'},'Location','best');
ylim([0 4.5])





