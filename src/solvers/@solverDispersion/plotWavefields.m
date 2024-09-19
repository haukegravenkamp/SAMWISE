function [figU,ax,wkMode] = plotWavefields(obj,msh,omegaK,opt,mat,bcd,linePlot)

% choose default options if not provided
if (nargin < 4) || ~isa(mat,'option')
    opt = option;
end
% remember whether materials are provided
% (required for acoustic displacement)
if (nargin < 5) || ~isa(mat,'material')
    materialProvided = false;
else
    materialProvided = true;
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if (nargin < 6) || ~isa(bcd,'bc')                                           % boundary conditions not provided (required for unbounded domains)
    bcd = bc;                                                               % set default
end

if (nargin < 7) || isempty(linePlot)                                        % flag for creating line plot instead of surface
    linePlot = false;                                                       % set default
end


%% interface
% solution
omegaPlot      = omegaK(1);
kPlot          = omegaK(2);
u              = obj.u;
omega          = obj.omega;
k              = obj.k;
pdes           = msh.pdes;
nPde           = msh.nPde;
sizes          = msh.size;
pde2elements   = msh.pde2elements;
materialNumber = msh.materialNumber;
eleOrder       = msh.eleOrder;
coordV         = msh.coordV;
dofsGlobal     = msh.dofsGlobal;
connectivity   = msh.connectivity;
nEle           = msh.nEle;
material2PDE   = msh.material2PDE;
pdeDofs        = msh.pdeDofs;

% options
plotFontSizes                = opt.plotting.plotFontSizes;
distortedGrid                = opt.plotting.distortedGrid;
acousticDisplacement         = opt.plotting.acousticDisplacement;
wavefieldPointsPerWavelength = opt.plotting.wavefieldPointsPerWavelength;
unboundedModel               = opt.model.unboundedModel;
waveFieldLx                  = opt.plotting.waveFieldLx;
waveFieldLz                  = opt.plotting.waveFieldLz;
waveFieldLineStyle           = opt.plotting.waveFieldLineStyle;
distortionFactor             = opt.plotting.distortionFactor;

% integration table
itg = inttable(msh.maxOrder);                                               % get integration table

thickness = sum(sizes);                                                     % total thickness of layered structure
lengthPlot = waveFieldLx*thickness;                                         % length in propagation direction to plot
lengthUnb = waveFieldLz*thickness;                                          % length in direction of unbounded domains coupled to surface

if isempty(u)
    figU = []; 
    ax = [];
    wkMode = [];
    warning('no eigenvectors provided for plotting wave fields.')
    return
end


%% find requested mode
[uPlot,wkMode,indwk] = getRequestedMode(omega,k,u,omegaPlot,kPlot);

%% rotate eigenvector arbitrarily to avoid purely real or imaginary components
uPlot = uPlot*exp(1i*pi/4);

%% check for coupling to unbounded domain
[coordV,connectivity,dofsGlobal,materialNumber,sizes,isUnbounded,pde2elements,unbBTa]=...
    includeUnboundedDomain(bcd,nEle,coordV,connectivity,...
    dofsGlobal,materialNumber,material2PDE,nPde,pde2elements,...
    sizes,lengthUnb,unboundedModel);

%% find smallest wavenumber to determine grid size
if linePlot                                                                 % standard graphs of mode shapes
    gridSize = min(msh.size./(msh.eleOrder*10));                             % choose grid size based on element order and thickness 
else
    gridSize = getGridSize(obj,wkMode(2),indwk,wavefieldPointsPerWavelength,unbBTa);
end


%% plot
figU = figure;
axis off
hold all

if linePlot
    nDivX = 1;
    lengthPlot = 0;
else
    nDivX = ceil(lengthPlot/gridSize);
end

rectID{nEle} = 0;                                                           % handles for plotting colored rectangles indicating element borders

xMin  = min(coordV(:,1));
xMax  = xMin + lengthPlot;
zMin  = min(coordV(:,3));
zMax  = max(coordV(:,3));
diffz = abs(zMax-zMin);

cols = cell(1,nPde);                                                        % initialize colormaps for each PDE

maxUpde = max(abs(uPlot));                                                  % maximum value of current EV for normalization
uPlot = uPlot/maxUpde;                                                      % normalize eigenvector

% normalize grid distortion
% default: use all dofs
maxDist = max(abs(uPlot));                                                  % maximum absolut displacement
normDist = distortionFactor*gridSize / maxDist;                             % absolute distortion factor
% in case of several PDEs, use only elastic dofs
for iPde = 1:nPde
    if isa(pdes(iPde),'pdeElasticity')
        elasDofs = nonnan(pdeDofs{iPde});
        maxDist = max(abs(uPlot(elasDofs(:))));
        normDist = distortionFactor*gridSize / maxDist;
    end
end
if linePlot                                                                 % line plot
    normDist = abs(zMax-zMin)/maxDist*0.3;                                  % normalize to fraction of overall thickness
end

% The normalization of grid distortion is currently not adequate if there
% are only acoustic media. As the acoustic displacement is computed later,
% we cannot predict what the maximum distortion will be. In those cases,
% the user should change the opt.plotting.distortionFactor to make the
% distortion visible.

for iPde = 1:nPde                                                           % loop PDEs

    pdeC = pdes(iPde);
    pdeName = class(pdeC);

    if ~linePlot
        ax(iPde) = axes;                                                        % create new axes object
        % each PDE can be plotted in a separate coordinate system, mainly to
        % use different colormaps to distinguish between elastic and acoustic
        % waves
    end
    hold all

    [cols,eleNos,dof,~,variableLabel] = ...
        generalPDEproperties(pdes,cols,iPde,pde2elements,...
        distortedGrid,acousticDisplacement,materialProvided);               % sort out some properties of the PDE

    maxDistZ = 0;                                                           % for keeping track of maximum vertical distortion
    for iE = 1:numel(eleNos)                                                % loop elements

        iEle = eleNos(iE);                                                  % global element number
        eleSize = sizes(iEle);                                              % current element size

        matEle = mat(materialNumber(iEle));                                 % materials

        eleDofs = nonnan(dofsGlobal(iEle,:));                               % dofs of current element
        vertexNumbers = connectivity(iEle,1:2);                             % vertex numbers (element end points)
        coordVele = coordV(vertexNumbers,:);                                % coordinates of vertices

        uEle = uPlot(eleDofs(:));                                           % part of eigenvector
        nDivZ = round(eleSize/gridSize);                                    % number of grid points in current element
        eta = linspace(-1,1,nDivZ);                                         % local coordinates

        [xyz,~,Jdet] = coordinateTransform(msh,2,coordVele,eta);              % global coordinates and Jacobian determinant
        
        % coordinateTransform(obj,nV,coordNodes,eta,N,Nd1,Nd2,~elementIsCurved(iEle))
        
        coordStart = coordVele(1,1,1);                                      % first point of element
        xEle1D = linspace(coordStart,coordStart+lengthPlot,nDivX);          % x values
        zEle1D = xyz(:,3);                                                  % y values
        [xEle,zEle] = meshgrid(xEle1D,zEle1D);                              % grid of x,y values

        % get solution at x- and y-values
        if isUnbounded(iEle)
            [uxz,pxz,sigxz,Pxz] = pointValuesUnbounded(obj,uEle,wkMode,indwk,xEle1D,zEle1D,nDivX,nDivZ,matEle,unbBTa(iEle),pdeName);
        else
            [uxz,pxz,sigxz,Pxz] = pointValuesElements(uEle,wkMode,Jdet,xEle1D,eta,eleOrder(iEle),itg,dof,nDivX,nDivZ,matEle,pdeName);
        end

        if distortedGrid                                                    % displacements used to distort grid
            xEle = xEle + uxz(:,:,1)*normDist;
            zEle = zEle + uxz(:,:,3)*normDist;
        end

        if ~acousticDisplacement && isa(pdeC,'pdeAcoustics')                % plot pressure, not displacement amplitude for acoustic problems
            amplitudePlot = vecnorm(pxz,2,3);                               % overwrite plot values by displacement amplitude, even if primary variable is pressure
        else
            amplitudePlot = vecnorm(uxz,2,3);                               % displacement amplitude for plotting
        end

        if linePlot
            ax(1)=subplot(1,3,1);
            [~,rectID{iEle,1}] = plotGraphs(uxz,cols{iPde},coordVele,zEle1D,maxDist,normDist,'u',(iPde == nPde)&&(iE == numel(eleNos)),opt);

            ax(2)=subplot(1,3,2);
            if ~isempty(sigxz)
                plotArray = sigxz;
                varName = '\sigma';
            else
                plotArray = pxz;
                varName = 'p';
            end
            [~,rectID{iEle,2}] = plotGraphs(plotArray,cols{iPde},coordVele,zEle1D,maxDist,normDist,varName,iE == numel(eleNos),opt);

            ax(3)=subplot(1,3,3);
            [~,rectID{iEle,3}] = plotGraphs(Pxz,cols{iPde},coordVele,zEle1D,maxDist,normDist,'P',(iPde == nPde)&&(iE == numel(eleNos)),opt);

        else
            surf(xEle,zEle,amplitudePlot,'FaceColor','interp','LineStyle',waveFieldLineStyle);
            hold all
            if ~isUnbounded(iEle) || unbBTa(iEle) == 2
                maxAmp = max(amplitudePlot(1,:));
                lineCol = [0 0 0];
                plot3(xEle(1,:),zEle(1,:),maxAmp*ones(size(xEle(1,:))),'-','Linewidth',1.5,'Color',lineCol)
            end
            if ~isUnbounded(iEle) || unbBTa(iEle) == 1
                maxAmp = max(amplitudePlot(end,:));
                lineCol = [0 0 0];
                plot3(xEle(end,:),zEle(end,:),maxAmp*ones(size(xEle(1,:))),'-','Linewidth',1.5,'Color',lineCol)
            end
            zlabel(pdes(iPde).variableName)
            cb = colorbar;
            cb.Label.String = variableLabel;
            axis equal;
        end
        if ismatrix(uxz)
            maxDistZ = max(maxDistZ,max(uxz(:,2)));
        elseif ndims(uxz)==3
            maxDistZ = max(maxDistZ,max(max(uxz(:,:,2))));
        end
    end

    grid off
    figureStandardSettings(plotFontSizes(1),plotFontSizes(2),1)

end

%% set axis properties
% especially when there is more than one axis and colorbar

if linePlot
    for iAx = 1:3
        subplot(1,3,iAx)
        ax(iAx).XLimitMethod = 'tight';
        % axis tight
        xL = ax(iAx).XLim;
        dxL = abs(diff(xL));
        for iR = 1:size(rectID,1)
            rectID{iR,iAx}.Position(1) = xL(1) - 0.5*dxL;
            rectID{iR,iAx}.Position(3) = 2*dxL;
        end
        ax(iAx).XLim=[xL(1) - 0.25*dxL, xL(2) + 0.25*dxL];
        ax(iAx).YLim = [zMin,zMax] + [-1,1]*diffz*0.05;
        grid off
        figureStandardSettings(plotFontSizes(1),plotFontSizes(2),1)
        axis off
    end
else
    ax(1).Colorbar.Location = 'northoutside';
    xL = [xMin,xMax] + [-1,1]*diffz*0.05;
    yL = [zMin,zMax] + [-1,1]*diffz*0.05 + [-1,1]*maxDistZ*normDist*distortedGrid;
    ax(1).XLim = xL;
    ax(1).YLim = yL;
    axPosition = ax(1).Position;
    wholeWidth = ax(1).Colorbar.Position(3);
    widthColorbar = wholeWidth/nPde*0.95;
    gapSize = (wholeWidth - nPde * widthColorbar)/(nPde+1);
    posXfirst = ax(1).Colorbar.Position(1);

    for iPde = 1:numel(ax)
        try
            ax(iPde).XLabel.String = '$x$';
            ax(iPde).YLabel.String = '$z$';
            ax(iPde).Colorbar.Location = 'northoutside';


            posXcolormap = widthColorbar*(iPde-1) + posXfirst + gapSize*iPde;
            ax(iPde).Colorbar.Position(1) = posXcolormap;
            ax(iPde).Colorbar.Position(3) = widthColorbar;
            ax(iPde).Colorbar.Position(2) = ax(1).Colorbar.Position(2);
            ax(iPde).Colorbar.Position(4) = ax(iPde).Colorbar.Position(4)*0.6;

            ax(iPde).XLim = xL;
            ax(iPde).YLim = yL;
            ax(iPde).Position = axPosition;

            ax(iPde).Color = 'none';
            ax(iPde).View = [0 90];
            ax(iPde).Colormap = cols{iPde};
            if iPde>1
                axis off
            end
        catch
            warning('unable to set axis properties correctly')
        end
    end


end
end


function [le, rectID] = plotGraphs(uxz,cols,coordVele,zEle1D,maxDist,normDist,varName,isLastEle,opt)

uxz = squeeze(uxz);
prop = {'Linewidth',1.5,'MarkerSize',5};                                    % plot properties
lineStyles = {'-','--','-.'};
nCol = size(cols,1);                                                        % number of colors in current colormap
ndofC = size(uxz,2);                                                        % number of dofs
if ndofC == 1
    labels = {['$',varName,'$']};
else
    labels = {...
        ['$',varName,'_x$'],...
        ['$',varName,'_y$']',...
        ['$',varName,'_z$'],...
        ['$',varName,'_{yz}$'],...
        ['$',varName,'_{xz}$'],...
        ['$',varName,'_{xy}$']};
end


markers = {'>','o','^','diamond','pentagram','hexagram'};
colSq = cols(round(nCol*0.95),:);
if mean(colSq)<0.9
    colSq = colSq/mean(colSq);
end
if max(colSq)>0.999
    colSq = colSq-(max(colSq)-1);
end
rectID = rectangle('Position',[1e-6, coordVele(1,3),...
    1e-6 coordVele(2,3)-coordVele(1,3)],...
    'FaceColor',colSq,'LineWidth',1,'EdgeColor','k',...
    'HandleVisibility','off');
legend('Location','best')
hold all
plot([0;0],coordVele(:,3),'--k','HandleVisibility','off')

le(ndofC) = matlab.graphics.chart.primitive.Line;

indM = round(size(uxz,1)/(ndofC+1)*(ndofC:-1:1));

maxU = max(abs(nonnan(uxz(:))));

for iDof = 1:ndofC                                              % loop dofs
    uC = uxz(:,iDof);
    if ~opt.plotting.plotZeroGraphs &&  ~any(abs(uC)>1e-6*maxU)
        continue
    end
    colC = cols( round(nCol*iDof/(ndofC+1.5)),:);           % pick color from colormap
    plot(uC*normDist,zEle1D,prop{:},...
        'Color',colC,...
        'DisplayName',labels{iDof},...
        'LineStyle',lineStyles{1},...
        'HandleVisibility','off');
    le(iDof)=plot(uC(indM(iDof))*normDist,zEle1D(indM(iDof)),...
        'Color',colC,'Linewidth',1, ...
        'DisplayName',labels{iDof}, 'MarkerFaceColor','w',...
        'LineStyle','none','Marker',markers{iDof});
    if ~isLastEle
        le(iDof).HandleVisibility = 'off';
    end
end


end
