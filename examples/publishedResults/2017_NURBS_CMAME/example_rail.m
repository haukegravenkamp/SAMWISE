clc
clear
close all

%% Example 5.3 Curved geometry
% Gravenkamp, Hauke, Sundararajan Natarajan, and Wolfgang Dornisch. “On the
% Use of NURBS-Based Discretizations in the Scaled Boundary Finite Element 
% Method for Wave Propagation Problems.” Computer Methods in Applied 
% Mechanics and Engineering 315 (2017): 867–80. 
% https://doi.org/10.1016/j.cma.2016.11.030.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

%% materials
materials = material;                                                       % create array of materials
materials.elastic.rho = 7.85;
materials.elastic.cs = 3.2;
materials.elastic.nu = 1/3;

%% geometry
geom = userMesh;                                                            % geometry defined by mesh file
geom.fileName = 'rail2.dat';                                                % mesh file name
geom.scaleGeometry = [1, 1]* 1/2.54*100;                                    % 
% The geometry is defined by a 'userMesh', i.e., a mesh file created by a
% different software. So far, only Ansys files are supported, others (gmsh,
% Abaqus) will be provided soon. Rather than specifying details of the
% geometry, userMesh just requires the name of the mesh file.
% As a little peculiarity, this mesh file happened to be created in inches 
% and also too small by a factor of 10. Hence, we use an option to scale
% the entire geometry by a given factor. scaleGeometry is a 1x2 vector in
% case we need different scaling factors in the y and z directions.

%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;

%% solver
sol = solverDispersion;
sol.fMin = 1e-5;
sol.fMax = 0.01;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
figure
plot(sol,'all',opt);

%% compare with disperse
ref = load('rail2017_cp.mat');
figure
axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
hold all
plot(ref.f/1e3,ref.cp,'.k')
legend([axH{1}.Children(1),axH{1}.Children(end)],{'samwise','CMAME 2017'},'Location','best');
xlim([sol.fMin,sol.fMax])
ylim([0 12])

