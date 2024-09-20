% applies some standard settings to current figure

% Select a figure, then run the function


function mystdfig(fontS,fontS2,figRatio)

%% default values
if (nargin<1)||isempty(fontS)
    fontS=12;                       % default font size
end
if (nargin<2)||isempty(fontS2)
    fontS2=fontS*14/12;
end
if (nargin<3)||isempty(figRatio)
    figRatio = 0;
end


%% figure properties
set(gcf,'windowstyle','normal')
set(gcf,'Color',[1 1 1])
scrsize=get(0,'screensize');

if length(findobj(gcf,'type','axes')) > 2
    increaseFactor = 1.4;
    posi = [scrsize(3)/10 scrsize(4)/10];
else
    increaseFactor = 1;
    posi = [scrsize(3)/3 scrsize(4)/3];
end

xSizeDefault = 610;
ySizeDefault = 377;

xSize = increaseFactor * xSizeDefault;

if figRatio==2                                                              % golden ratio
    ySize = round(xSize/1.6180);
    axis normal
elseif figRatio==1                                                          % square
    ySize = ySizeDefault;   
    xSize = ySizeDefault;
    axis square
else                                                                        % default
    ySize = round(xSize * 420/610);
    axis normal
end

set(gcf,'position',[posi(1) posi(2) xSize ySize])

clear scrsize

set(gcf,'defaulttextinterpreter','latex')
set(gcf,'DefaultTextFontname', 'CMU Serif')
set(groot, 'DefaultLegendInterpreter', 'latex')


%% axis properties
ax=gca;

ax.Title.FontSize=fontS2;
ax.FontName='CMU Serif';
ax.FontSize = fontS;
ax.GridLineWidth = 0.5;
ax.TickLabelInterpreter = 'latex';
ax.XLabel.FontSize = fontS2;
ax.YLabel.FontSize = fontS2;
ax.ZLabel.FontSize = fontS2;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';
set(gca,'Linewidth',0.5)

box on

% colorbar
if ~isempty(ax.Colorbar)
    ax.Colorbar.Label.Interpreter = 'latex';
    ax.Colorbar.Title.Interpreter = 'latex';
    ax.Colorbar.Label.FontSize = fontS2;
    ax.Colorbar.Title.FontSize = fontS2;
    ax.Colorbar.FontName = 'CMU Serif';
    ax.Colorbar.TickLabelInterpreter='latex';
end

% legend
hLegend = findobj(gcf, 'Type', 'Legend');
if ~isempty(hLegend)
    for i=1:numel(hLegend)
        hLegend(i).FontName = 'CMU Serif';
        hLegend(i).Interpreter = 'latex';
        hLegend(i).FontSize = fontS;
    end
end

end
