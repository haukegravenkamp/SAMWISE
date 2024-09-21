
clc
clear
close all

%% example 6.1: Brass plate immersed in water
% Gravenkamp, Hauke, Bor Plestenjak, Daniel A Kiefer, and Elias Jarlebring.
% “Computation of Leaky Waves in Layered Structures Coupled to Unbounded
% Media by Exploiting Multiparameter Eigenvalue Problems.” Journal of Sound
% And Vibration, 2024. https://doi.org/10.1016/j.jsv.2024.118716.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% materials

materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % name for reading from database
materials(2).name = 'water';

%% geometry

geom = plate;                                                               % create geometry of type plate                                                             % layer number
geom.layers.thickness = 1;                                                  % thickness
geom.layers.material = 1;                                                   % material number

%% boundary conditions

bCond(1) = bcUnbounded;                                                     % coupling to unbounded domain
bCond(1).location = 'bottom';                                               % lower surface
bCond(1).material = 2;                                                      % material number of unbounded domain (refers to materials array above)

bCond(2) = bcUnbounded;                                                     % same for top surface
bCond(2).location = 'top';
bCond(2).material = 2;
% We create two boundary conditions of type bcUnbounded, one for the bottom
% and one for the top surface. In this case, both surfaces are coupled to
% the same infinite material (water). For other cases, including different
% materials and coupling to solids, see the other examples in this folder.

%% settings

opt = option;
opt.model.unboundedModel = 'exact';                                         % exact boundary conditions -> multipareig solver
% opt.model.unboundedModel = 'dashpot';                                     % approximate dashpot boundary conditions -> standard solver
% when the model involves coupling to an unbounded domain, we can either
% choose a fast but not necessarily accurate 'dashpot' approximation, or a
% rather expensive but exact formulation of the boundary condition.

opt.plotting.maxAttenuation = 2000;                                         % maximum attenuation of modes to plot
opt.plotting.waveFieldLx = 4;                                               % how long of a plate section to plot when displaying wave fields
opt.plotting.waveFieldLz = 1;                                               % extent of infinite domain in plot

opt.numerics.unboundedRemoveIncoming = ~true;                               % whether to filter out incoming waves in postprocessing
opt.numerics.removeNegativeAttenuation = ~true;                             % whether to filter out modes with negative attenuation
opt.numerics.unboundedRemoveThreshold = -1e-2;                              % relative threshold for finding incoming waves

%% solver
sol = solverDispersion;
sol.fMin = 0.001;
sol.fMax = 4;
sol.nSteps = 300;
% note that it's sufficient to initialize the standard solver. The subclass
% will be selection based on the boundary conditions and the unboundedModel
% option.

%% call program
[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(materials, geom, opt, bCond, sol);

%% plot wave fields (2d surfaces)

omKplot1 = [1*2*pi,0.25*2*pi];                                              % pick (k,omega) of mode to plot
opt.plotting.wavefieldPointsPerWavelength = 25;                             % number of plotting points per wave length
[~,~,omegaKmode1] = plotWavefields(sol,msh,omKplot1,opt,mat,bcd);           % plot requested mode
set(gcf,'Position',[647   357   840   519])
drawnow

omKplot2 = [3.5*2*pi,8];                                                    % same for second mode
opt.plotting.wavefieldPointsPerWavelength = 12;                             
[~,~,omegaKmode2] = plotWavefields(sol,msh,omKplot2,opt,mat,bcd);
set(gcf,'Position',[647   357   840   519])

%% plot mode shape (line plots)
plotModeshape(sol,msh,omKplot1,opt,mat,bcd);
plotModeshape(sol,msh,omKplot2,opt,mat,bcd);
drawnow

%% compare with stored solution, obtained using Kiefer's linearization

ref = load('plate_immersed_BrassWater.mat');                                % load reference solution

% define legend entry
if strcmp(opt.model.unboundedModel, 'dashpot')
    displName = 'dashpot';
elseif strcmp(opt.model.unboundedModel, 'exact')
    displName = 'MultiParEig';
end

plotParam = {'LineWidth',0.5,'LineStyle','none','Color',[0.17,0.49,0.63],'Marker','o','MarkerSize',6,'MarkerFaceColor',[1 1 1]*0.9,'DisplayName',displName};

figure
hold all
axH = plot(sol,'cp',opt,plotParam);
plot(ref.fcp,ref.cp,'.k','MarkerSize',5)
legend([axH{1}.Children(end),axH{1}.Children(1)],{displName,'linearization'},'Location','best');
ylim([0 10])
plotParamPoints = {'Color',myColor(2),'LineWidth',1.5,'MarkerSize',8,'HandleVisibility','off'};
fPl1 = real(omegaKmode1(1))/2/pi;
fPl2 = real(omegaKmode2(1))/2/pi;
plot(fPl1,real(omegaKmode1(1))/real(omegaKmode1(2)),plotParamPoints{:},'Marker','square');
plot(fPl2,real(omegaKmode2(1))/real(omegaKmode2(2)),plotParamPoints{:},'Marker','diamond');

figure
hold all
axH = plot(sol,'att',opt,plotParam);
plot(ref.fatt,ref.att,'.k','MarkerSize',5)
legend([axH{1}.Children(end),axH{1}.Children(1)],{displName,'linearization'},'Location','best');
ylim([0 2000])
plotParamPoints = {'Color',myColor(2),'LineWidth',1.5,'MarkerSize',8,'HandleVisibility','off'};
plot(fPl1,imag(omegaKmode1(2))*20/log(10)*1000,plotParamPoints{:},'Marker','square');
plot(fPl2,imag(omegaKmode2(2))*20/log(10)*1000,plotParamPoints{:},'Marker','diamond');

