function [axisHandles, plotHandles] = plotSolver(obj,plotVar,opt,plotParameters)

if (nargin < 2) || isempty(plotVar)
    plotVar = 'all';
end

if (nargin < 3) || isempty(opt)
    opt = option;
end
if (nargin < 4) || isempty(plotParameters)
    plotParameters = {...
        'Color',[0.17,0.49,0.63],...
        'Marker','.',...
        'MarkerSize',8,...
        'MarkerFaceColor','w',...
        'LineStyle','none'...
        };
end

%% interface
% solution
f     = obj.f;
omega = obj.omega;
k     = obj.k;
cp    = obj.cp;
cg    = obj.cg;
att   = obj.att;
% options
useAngularFrequency       = opt.plotting.angularFrequency;
omegaVSk                  = opt.plotting.omegaVSk;
labelWavenumber           = opt.plotting.labelWavenumber;
labelPhaseVelocity        = opt.plotting.labelPhaseVelocity;
labelGroupVelocity        = opt.plotting.labelGroupVelocity;
labelFrequency            = opt.plotting.labelFrequency;
labelCircularFrequency    = opt.plotting.labelCircularFrequency;
labelAttenuation          = opt.plotting.labelAttenuation;
plotFontSizes             = opt.plotting.plotFontSizes;
propagatingModesThreshold = opt.plotting.propagatingModesThreshold;
maxAttenuation            = opt.plotting.maxAttenuation;
omitNegativeWavenumber    = opt.plotting.omitNegativeWavenumber;
omitNegativeAttenuation   = opt.plotting.omitNegativeAttenuation;

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% check which results to plot
plotInd = false(4,1);                                                       % flag for plotting k, cp, cg, att
varNamesShort = {'k','cp','cg','att'};                                      % recognized variable names, abbreviation
varNames = {'wavenumber','phase velocity','group velocity','attenuation'};  % recognized variable names, full name
plotAll = false;
for i = 1:4                                                                 % loop variable names
    if strcmpi(plotVar,varNamesShort{i}) || strcmpi(plotVar,varNames{i})    % option is in list of variables
        plotInd(i) = true;                                                  % set flag
        break
    end
end
if ~any(plotInd)
    plotInd(:) = true;                                                      % plot everything if no specific flag was found
    plotAll = true;                                                         % flag for plotting everything
end


%% filter requested modes
if maxAttenuation == 0
    maxAttenuation = propagatingModesThreshold;
end
indR = (abs(att) > maxAttenuation);                                         % indices of modes with too large attenuation
if omitNegativeAttenuation
    indR = indR | (att < -maxAttenuation*1e-3);
end
if omitNegativeWavenumber
    indR = indR | (real(k) < -eps);
end
ind = ~indR;                                                                % indices of modes to keep

if numel(find(ind))/numel(ind) < 1/200
    warning('There are not many modes to plot. You may want to increase the maxAttenuation value.')
end
k(~ind) = nan;
cp(~ind) = nan;
if ~isempty(cg)
cg(~ind) = nan;
end
att(~ind) = nan;

for iw = 1:size(k,2)
    [~,ind] = sort(real(k(:,iw)));
    k(:,iw) = k(ind,iw);
    [~,ind] = sort(real(cp(:,iw)));
    cp(:,iw) = cp(ind,iw);
    if ~isempty(cg)
    [~,ind] = sort(real(cg(:,iw)));
    cg(:,iw) = cg(ind,iw);
    end
    [~,ind] = sort(real(att(:,iw)));
    att(:,iw) = att(ind,iw);
end
ind = ~all(isnan(k),2);
k = k(ind,:);
cp = cp(ind,:);
if ~isempty(cg)
cg = cg(ind,:);
end
att = att(ind,:);

kReal = real(k);
plotVars = {k,cp,cg,att};                                                   % put in cell

%% some settings for axis limits
% determine limits for phase velocity based on largest wave number
%
[~,indMax] = max(kReal(:,end),[],'all');
cpMax = omega(end)/kReal(indMax,end);
yLimCp = [0,3*cpMax];
yLimCg = yLimCp;
if yLimCp(2) < yLimCp(2) + eps
    yLimCp(2) = yLimCp(2) + 2*eps;
    yLimCg(2) = (1+eps)*yLimCg(2) + 2*eps;
end
if maxAttenuation == 0
    yLimAtt = [0, 1e-6];
else
    yLimAtt = [0, maxAttenuation];
end


%% plot
if useAngularFrequency
    fPlot = omega;
    fLabel = labelCircularFrequency;
else
    fPlot = f;
    fLabel = labelFrequency;
end
plotLabels = {labelWavenumber,labelPhaseVelocity,labelGroupVelocity,labelAttenuation};
ylimits = {'auto',yLimCp,yLimCg,yLimAtt};

plotHandles = {};
axisHandles = {};
for iPlot = 1:4                                                             % loop variables
    if ~plotInd(iPlot)                                                      % current variable not to be plotted
        continue                                                            % do nothing
    end
    if plotAll                                                              % plot all variables
        axisHandles{iPlot} = subplot(2,2,iPlot);                            % subplot axis
    else                                                                    % only one variable
        axisHandles{1} = gca;                                               % only one set of axes
    end
    xPlot = fPlot;                                                          % frequency on x-axis usually
    yPlot = plotVars{iPlot}';                                                % current variable on y-axis
    xlab = fLabel;                                                          % frequency label
    ylab = plotLabels{iPlot};                                               % variable label
    if (iPlot == 1) && omegaVSk                                             % wavenumber plot and omega vs k requested
        xPlot = yPlot;                                                      % switch axes
        yPlot = repmat(fPlot,size(xPlot,2),1);
        xlab = ylab;                                                        % switch labels
        ylab = fLabel;
    else
        xPlot = repmat(xPlot,size(yPlot,2),1);
    end

    plotHandles{iPlot} = plotVariable(xPlot,yPlot.',xlab,ylab,ylimits{iPlot},plotParameters,plotFontSizes);

end

    function plHandle = plotVariable(x,y,xlab,ylab,yl,plotParameters,fontS)
        plHandle  = plot(real(x),real(y),plotParameters{:});
        xlabel(xlab);
        ylabel(ylab);
        figureStandardSettings(fontS(1),fontS(2),true)
        axis tight
        ylim(yl)
    end




end
