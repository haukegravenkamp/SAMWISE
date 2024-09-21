%% example of a brass plate
% computes dispersion curves of a 1mm-thick brass plate with traction-free
% surfaces and compares the results with those obtained using 'disperse'

%% materials
mat = material;                                                             % create array of materials
mat.name = 'brass';                                                         % name for reading from database
% Here, we only need one material. If we provide a name, it will be checked
% whether this name exists in the database
% \src\materialDefinition\materialDatabase. If so, the properties stored
% there will be used, any other properties will be ignored. If the name is
% not found in the database, the code will assume the user has provided
% properties here; if those are incomplete, we will get an error. In this
% example, we are happy with the parameters stored in the database under
% the name 'brass'.

%% geometry
geom = plate;                                                               % create geometry of type plate
% see example_minimal

geom.layers.thickness = 1;                                                  % assign thickness of 1 mm
geom.layers.material = 1;                                                   % material number
% This is the index in the array 'mat' created above; here, we only have
% one material.

%% settings 
opt = option;                                                               % create options
% option is a class used for storing various settings. It is initialized
% with default values. The options are divided into the categories 
% 'model', 'numerics', 'postprocessing', 'plotting'. To find out the
% meaning of the options, look into the corresponding files in \src\options
opt.plotting.omitNegativeWavenumber = true;
% this option prevents solutions with negative wavenumber being plotted

%% solver
sol = solverDispersion;                                                     % create solver
% see example_minimal. But here, we explicitly choose the minimum and
% maximum frequency and number of steps
sol.fMin = 0;
sol.fMax = 7;
sol.nSteps = 100;

%% call program
[geo, mat, bcd, sol, opt, res, msh] = samwise(mat, geom, opt, sol);
% Here, we pass four opjects to the program, again in arbitrary order

%% plot
plot(sol,'all',opt);
% plot dispersion curves
% Note that we pass the options to the plot function so that it considers
% the omitNegativeWavenumber setting. The keyword 'all' means that we want
% to plot wavenumber, phase velocity, group velocity, and attenuation. 
% We can choose them individually by replacing 'all' by any of 
% 'k', 'cp', 'cg', 'att', e.g., plot(sol,'cg',opt);

%% compare with disperse
ref = load('plateBrass_disperse.mat');                                      % load previously stored solution obtained by disperse software
figure
hold all
plot(ref.f,ref.cp,'-k','LineWidth',1.5)                                     % plot disperse solution

axH = plot(sol,'cp',opt,{'Markersize',3,'Color',[0.1725,0.4902,0.6275],'LineStyle','none','Marker','o','MarkerFaceColor','w'});
legend([axH{1}.Children(1),axH{1}.Children(end)],{'samwise','disperse'},'Location','best');
% plot samwise solution
% Here, we set the flag 'cp' to only plot the phase velocity. We also
% provide a bunch of optional plot properties, which are the same as for
% Matlab's standard plot function.
% Out plot function returns a cell of handles to each axis, which is used
% here to modify the legend.

