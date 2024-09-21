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
geom = cylinder;                                                            % initialize cylinder geometry
geom.layers.thickness = 6;                                                  % thickness
geom.layers.material  = 1;                                                  % material number
geom.r0 = 3;
% Here, the geometry is a cylinder, defined by the wall thickness (!),
% inner radius, and material number

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 1.9;
sol.nSteps = 100;

%% options
opt = option;
opt.model.circumferentialOrders = 0;
% For a cylinder, we need to define the order n of modes in the
% circumferential direction that we wish to compute. This refers to the
% modes proportional to exp(i*n*theta) where theta is the angle in a
% cylindrical coordinate system.
% The default is n = 0:1
% Note that n=0 includes both longitudinal and torsional modes. The
% reference solution only shows longitudinal modes.

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(mat, geom, sol,opt);

%% compare with Global Matrix Method
ref = load('cylinderPPN.mat');
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1)
axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
legend([axH{1}.Children(1),axH{1}.Children(end)],{'samwise','GMM'},'Location','best');
ylim([0 4.5])





