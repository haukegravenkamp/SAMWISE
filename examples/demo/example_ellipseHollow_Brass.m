clc
clear
close all

%% materials
materials = material;                                                       % create array of materials
materials(1).name = 'brass';                                                % name for reading from database

%% geometry
geom = ellipseHollow;                                                            % full cylinder
geom.Ao = 4;                                                                 % outer radius
geom.Ai = 3;                                                                 % inner radius
geom.Bo = 2;                                                                 % outer radius
geom.Bi = 1;                                                                 % inner radius


%% settings 
opt = option;
opt.plotting.omitNegativeWavenumber = true;
opt.numerics.eleOrderFunction = @(w)ceil(3+0.5*w);
opt.numerics.nEigenvaluesDefault = 50;

%% solver
sol = solverDispersion;
sol.fMin = 1e-3;
sol.fMax = 0.6;
sol.nSteps = 50;


%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(materials, geom, opt, sol);

%% plot
figure
plot(sol,'all',opt);

