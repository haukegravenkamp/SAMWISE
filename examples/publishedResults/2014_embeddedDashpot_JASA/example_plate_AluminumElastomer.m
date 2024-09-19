
%% Brass plate coupled on one side to unbounded Teflon

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
clc
clear 
close all

%% Example III.A Plate attached to infinite medium
% Gravenkamp, Hauke, Carolin Birk, and Chongmin Song. “Computation of 
% Dispersion Curves for Embedded Waveguides Using a Dashpot Boundary 
% Condition.” The Journal of the Acoustical Society of America 135, 
% no. 3 (2014): 1127–38. https://doi.org/10.1121/1.4864303.


%% materials
% according to Castaings et al. 'Finite element model for waves guided 
% along solid systems of arbitrary section coupled to infinite solid
% media', JASA 2008
%
% parameters in the paper are given by C11 and C66 = G, convert to G, nu
C11 = 112;
C66 = 27;
rho = 2.78;
nu = (2*C66-C11)/(2*C66-2*C11);

mat = material;                                                             % create array of materials
mat(1).name = 'aluminum';                                                   % name 
mat(1).elastic.G   = C66;                                                    % shear modulus
mat(1).elastic.rho = rho;                                                  % density
mat(1).elastic.nu  = nu;                                                % Poisson's ratio

C11 = 1+0.1*1i;
C66 = 0.3*(1+0.1i);
rho = 1;
nu = (2*C66-C11)/(2*C66-2*C11);

mat(2).name = 'elastomer';           
mat(2).elastic.G   = C66;
mat(2).elastic.rho = rho;
mat(2).elastic.nu  = nu;

%% geometry

geom = plate;                                            % create geometry of type plate
geom.assumption2D = 'none';                                                 % all three dofs
geom.layers.thickness = 4;                                                  % thickness
geom.layers.material  = 1;                                                  % material number of the layer


%% boundary conditions

bCond = bcUnbounded;                                                        % coupling to unbounded domain
bCond.location = 'top';                                                     % upper surface
bCond.material = 2;                                                         % material number of unbounded domain (refers to materials array above)

%% settings 

opt = option;
% opt.model.unboundedModel = 'exact';                                       % exact boundary conditions -> multipareig solver
opt.model.unboundedModel = 'dashpot';                                       % approximate dashpot boundary conditions -> standard solver
opt.plotting.maxAttenuation = 200;                                          % maximum attenuation, just for filtering/plotting

%% solver

sol = solverDispersion;
sol.fMin = 0.0;
sol.fMax = 0.4;
sol.nSteps = 100;

%% call program

[geo, mat, bcd, sol, opt, res, msh, glb] = samwise(mat, geom, opt, bCond, sol);

%% compare with absorbing region
ref = load('plate_AluminumElastomer_Castaings.mat');
figure
hold all
axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
plot(ref.fcp,ref.cp,'xk','LineWidth',1)
legend([axH{1}.Children(end),axH{1}.Children(1)],{'samwise','absorbing region'},'Location','best');
ylim([0 6])

figure
hold all
axH = plot(sol,'att',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
plot(ref.fatt,ref.att,'xk','LineWidth',1)
legend([axH{1}.Children(end),axH{1}.Children(1)],{'samwise','absorbing region'},'Location','best');



