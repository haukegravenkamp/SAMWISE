clc
clear 
close all

%% Example 5.1 Functionally graded material
% Gravenkamp, Hauke, Chongmin Song, and Jens Prager. “A Numerical Approach 
% for the Computation of Dispersion Relations for Plate Structures Using 
% the Scaled Boundary Finite Element Method.” Journal of Sound and 
% Vibration 331 (2012): 2543–57. https://doi.org/10.1016/j.jsv.2012.01.029.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials
rho1 = 7.19;
rho2 = 3.9;
lambda1 =  74.2;
lambda2 = 138;
G1 = 102.5;
G2 = 118.11;
p  = 10;
materials = material;                                                       % create array of materials
materials.name = 'FGM';                                                     % arbitrary name
materials.elastic.rho    = @(y,z) rho1*(1-z.^p) + rho2*(z.^p);              % define materials as function handles
materials.elastic.G      = @(y,z) G1*(1-z.^p) + G2*(z.^p);
materials.elastic.lambda = @(y,z) lambda1*(1-z.^p) + lambda2*(z.^p);
% note: when the material properties are defined as functions, the wave
% velocities are not uniquely defined for each material, and the element
% order is not chosen automatically. Hence, it will use the default element
% order given by opt.numerics.eleOrderDefault, see below.

%% geometry
geom = plate;                                                               % initialize plate geometry
geom.layers.thickness = 1;                                                  % layer thickness
geom.layers.material  = 1;                                                  % FGM materil defined above

%% solver
sol = solverDispersion;
sol.fMin   = 0;
sol.fMax   = 15;
sol.nSteps = 100;

%% options
opt = option;
opt.plotting.angularFrequency = ~true;
opt.numerics.eleOrderDefault = 15;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, sol);

%% compare with reference solution
figure
hold all
ref = load('FGM_Cao.mat');                                                  % load reference solution
plot(ref.k,ref.cp,'-k','LineWidth',1)
figureStandardSettings

% the plot cp vs. k from the literature is not standard, so we do it
% manually
thre = opt.plotting.propagatingModesThreshold;
indComplex = (abs(sol.att)>(thre));
kPlot = real(sol.k);
cPlot = real(sol.cp);
kPlot(indComplex) = nan;
cPlot(indComplex) = nan;

plot(kPlot,cPlot,'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w')
a = gca;
legend([a.Children(1),a.Children(end)],{'samwise','Cao et al. 2010'},'Location','best');
xlim([0 10])
ylim([0 10])
xlabel('$k$')
ylabel('$c_p$')
