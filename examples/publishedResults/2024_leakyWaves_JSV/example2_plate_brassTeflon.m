
clc
clear
close all

%% example 6.2: Brass plate coupled to infinite Teflon
% Gravenkamp, Hauke, Bor Plestenjak, Daniel A Kiefer, and Elias Jarlebring.
% “Computation of Leaky Waves in Layered Structures Coupled to Unbounded
% Media by Exploiting Multiparameter Eigenvalue Problems.” Journal of Sound
% And Vibration, 2024. https://doi.org/10.1016/j.jsv.2024.118716.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials

materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % name for reading from database
materials(2).name = 'teflon';

%% geometry

geom = plate;                                                               % create geometry of type plate
geom.assumption2D = 'none';                                                 % all three dofs
geom.layers.thickness = 1;                                                  % thickness
geom.layers.material  = 1;                                                  % material number


%% boundary conditions

bCond = bcUnbounded;                                                        % coupling to unbounded domain
bCond.location = 'top';                                                     % upper surface
bCond.material = 2;                                                         % material number of unbounded domain (refers to materials array above)

%% settings

opt = option;
opt.model.unboundedModel = 'exact';                                         % exact boundary conditions -> multipareig solver
% opt.model.unboundedModel = 'dashpot';                                     % approximate dashpot boundary conditions -> standard solver

opt.plotting.maxAttenuation = 4000;
opt.numerics.unboundedRemoveIncoming = true;
opt.numerics.removeNegativeAttenuation = true;
opt.numerics.unboundedUsePoyntingVec = ~true;
opt.numerics.unboundedRemoveThreshold = 1e-3;

%% solver

sol = solverDispersion;
sol.fMin = 0.01;
sol.fMax = 7;
sol.nSteps = 281;

%% call program

[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(materials, geom, opt, bCond, sol);

%% plot wave fields


omKplot1 = [1*2*pi,1.6];
omKplot2 = [3*2*pi,3*2*pi/7.34787];

opt.plotting.wavefieldPointsPerWavelength = 10;
opt.plotting.waveFieldLx = 3;
opt.plotting.waveFieldLz = 2;
[~,~,omegaKmode1] = plotWavefields(sol,msh,omKplot1,opt,mat,bcd);
set(gcf,'Position',[647   357   840   519])
colormap bone
drawnow

opt.plotting.wavefieldPointsPerWavelength = 5;
opt.plotting.waveFieldLx = 3;
opt.plotting.waveFieldLz = 2;
[~,~,omegaKmode2] = plotWavefields(sol,msh,omKplot2,opt,mat,bcd);
set(gcf,'Position',[647   357   840   519])
colormap bone
drawnow

%% plot mode shape (line plots)
opt.plotting.wavefieldPointsPerWavelength = 50;
plotModeshape(sol,msh,omKplot1,opt,mat,bcd);
plotModeshape(sol,msh,omKplot2,opt,mat,bcd);
drawnow

%% compare with stored solution, obtained using disperse

ref = load('plate_embedded_BrassTeflon.mat');                                % load reference solution

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
ylim([0 4000])
plotParamPoints = {'Color',myColor(2),'LineWidth',1.5,'MarkerSize',8,'HandleVisibility','off'};
plot(fPl1,imag(omegaKmode1(2))*20/log(10)*1000,plotParamPoints{:},'Marker','square');
plot(fPl2,imag(omegaKmode2(2))*20/log(10)*1000,plotParamPoints{:},'Marker','diamond');
