
clc
clear
close all

%% example 6.3: Layered plate between a fluid and an elastic half-space
% Gravenkamp, Hauke, Bor Plestenjak, Daniel A Kiefer, and Elias Jarlebring.
% “Computation of Leaky Waves in Layered Structures Coupled to Unbounded
% Media by Exploiting Multiparameter Eigenvalue Problems.” Journal of Sound
% And Vibration, 2024. https://doi.org/10.1016/j.jsv.2024.118716.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials

materials = material;                                                       % create array of materials

materials(1).name = 'brass';                                                % name for reading from database
materials(2).name = 'titanium';
materials(3).name = 'teflon';
materials(4).name = 'motorOil';


%% geometry

geom = plate;                                                               % create geometry of type plate
geom.assumption2D = 'planeStrain';

geom.layers(1).thickness = 1;                                              % thickness
geom.layers(1).material = 2;                                               % material number

geom.layers(2).thickness = 1;                                              % thickness
geom.layers(2).material = 1;                                               % material number

geom.layers(3).thickness = 1;                                              % thickness
geom.layers(3).material = 2;                                               % material number


%% boundary conditions

bCond(1) = bcUnbounded;                                                     % coupling to unbounded domain
bCond(1).location = 'bottom';                                               % upper surface
bCond(1).material = 3;                                                      % material number of unbounded domain (refers to materials array above)

bCond(2) = bcUnbounded;                                                     % coupling to unbounded domain
bCond(2).location = 'top';                                                  % upper surface
bCond(2).material = 4;                                                      % material number of unbounded domain (refers to materials array above)

%% settings

opt = option;
opt.model.unboundedModel = 'exact';                                         % exact boundary conditions -> multipareig solver
% opt.model.unboundedModel = 'dashpot';                                     % approximate dashpot boundary conditions -> standard solver

opt.plotting.maxAttenuation = 3000;

opt.numerics.unboundedRemoveIncoming = true;
opt.numerics.removeNegativeAttenuation = true;
opt.numerics.unboundedUsePoyntingVec = ~true;
opt.numerics.unboundedRemoveThreshold = -1e-2*0;

%% solver

sol = solverDispersion;
sol.fMin = 0.000;
sol.fMax = 3;
sol.nSteps = 121;

%% call program

[geom, materials, bCond, sol, opt, res, msh, glb] = samwise(materials, geom, opt, bCond, sol);


%% plot wave fields
omPlot1 = 2*pi;
kPlot1 = omPlot1/2.92 + 1i*att2k(660);
omKplot1 = [omPlot1,kPlot1];

omPlot2 = 2.5*2*pi;
kPlot2 = omPlot2/3.94 + 1i*att2k(853);
omKplot2 = [omPlot2,kPlot2];

opt.plotting.wavefieldPointsPerWavelength = 7;
opt.plotting.waveFieldLx = 2;
opt.plotting.waveFieldLz = 2/3;
opt.plotting.distortionFactor = 1;
figure
[~,~,omegaKmode1] = plotWavefields(sol,msh,omKplot1,opt,materials,bCond);
set(gcf,'Position',[647   357   840   519])

drawnow
opt.plotting.wavefieldPointsPerWavelength = 5;
opt.plotting.waveFieldLx = 2;
opt.plotting.waveFieldLz = 2/3;
opt.plotting.distortionFactor = 1.5;
[~,~,omegaKmode2] = plotWavefields(sol,msh,omKplot2,opt,materials,bCond);
set(gcf,'Position',[647   357   840   519])
drawnow

%% plot mode shape (line plots)
plotModeshape(sol,msh,omKplot1,opt,materials,bCond);
plotModeshape(sol,msh,omKplot2,opt,materials,bCond);
drawnow

%% compare with stored solution, obtained using disperse

ref = load('plate_compositeTeflonOil.mat');                                % load reference solution

% define legend entry
if strcmp(opt.model.unboundedModel, 'dashpot')
    displName = 'dashpot';
elseif strcmp(opt.model.unboundedModel, 'exact')
    displName = 'MultiParEig';
end

plotParam = {'LineStyle','none','Color',[0.17,0.49,0.63],'Marker','o','MarkerSize',5,'MarkerFaceColor',[1 1 1]*0.85,'DisplayName',displName};

figure
hold all
axH = plot(sol,'cp',opt,plotParam);
plot(ref.fcp,ref.cp,'-k','LineWidth',1.5)
legend([axH{1}.Children(end),axH{1}.Children(1)],{displName,'disperse'},'Location','best');
ylim([0 10])
plotParamPoints = {'Color',myColor(2),'LineWidth',1.5,'MarkerSize',8,'HandleVisibility','off'};
fPl1 = real(omegaKmode1(1))/2/pi;
fPl2 = real(omegaKmode2(1))/2/pi;
plot(fPl1,real(omegaKmode1(1))/real(omegaKmode1(2)),plotParamPoints{:},'Marker','square');
plot(fPl2,real(omegaKmode2(1))/real(omegaKmode2(2)),plotParamPoints{:},'Marker','diamond');

figure
hold all
axH = plot(sol,'att',opt,plotParam);
plot(ref.fatt,ref.att,'-k','LineWidth',1.5)
legend([axH{1}.Children(end),axH{1}.Children(1)],{displName,'disperse'},'Location','best');
ylim([0 2000])
plotParamPoints = {'Color',myColor(2),'LineWidth',1.5,'MarkerSize',8,'HandleVisibility','off'};
plot(fPl1,imag(omegaKmode1(2))*20/log(10)*1000,plotParamPoints{:},'Marker','square');
plot(fPl2,imag(omegaKmode2(2))*20/log(10)*1000,plotParamPoints{:},'Marker','diamond');
