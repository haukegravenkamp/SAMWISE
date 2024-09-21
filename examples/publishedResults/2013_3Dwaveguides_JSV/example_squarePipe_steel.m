clc
clear
close all

%% Example 6.2 Square pipe
% Gravenkamp, Hauke, Hou Man, Chongmin Song, and Jens Prager. “The 
% Computation of Dispersion Relations for Three-Dimensional Elastic 
% Waveguides Using the Scaled Boundary Finite Element Method.” Journal of 
% Sound and Vibration 332 (2013): 3756–71. 
% https://doi.org/10.1016/j.jsv.2013.02.007.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

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
% This geometry is an 'L-shape' - in other words, a quarter of a
% rectangular pipe. It is defined by the width Ly, height Lz, and the wall
% thicknesses hy, hz. Open the file Lshape.m for clarification.

%% boundary conditions
bcd(1) = bcFixed;
bcd(1).location = 'bottom';
bcd(1).directions = 3;                                                      % vertical direction (symmetric)
bcd(2) = bcFixed;
bcd(2).location = 'left';
bcd(2).directions = 2;                                                      % horizontal direction (symmetric)
% We define two boundary conditions, both of them fixing some degrees of
% freedom. The first one fixes all z-components at the bottom of the
% cross-section, the secton one fixes the y-direction at the left. This
% corresponds to symmetric boundary conditions on the y- and z-axes.

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;

%% solver
sol = solverDispersion;
sol.fMin = 0;
sol.fMax = 0.1;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(materials, geom, opt, sol, bcd);

%% plot
figure
plot(sol,'all',opt);

%% compare with old results
ref = load('squarePipe_Sorohan.mat');
figure
hold all
axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
plot(ref.f/1000,ref.cp,'.k','MarkerSize',12)
legend([axH{1}.Children(end),axH{1}.Children(1)],{'samwise','Sorohan et al. 2010'},'Location','best');
ylim([0 10])
xlim([0,0.1])
