
function figureStandardSettings(fontS,fontS2,figRatio)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com

if nargin<1
    fontS = 12;                       % default font size
end
if nargin<2
    fontS2 = 14;
end
if nargin<3
    figRatio = 1;
end


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

% xSizeDefault = 600;
xSizeDefault = round(scrsize(3)/2);
xSize = increaseFactor * xSizeDefault;

 
if figRatio
    ySize = round(xSize/1.6180);
else
    ySize = round(xSize * 420/610);
end
    
set(gcf,'position',[posi(1) posi(2) xSize ySize])


clear scrsize

set(gcf,'defaulttextinterpreter','latex')
set(gcf,'DefaultTextFontname', 'CMU Serif')
set(groot, 'DefaultLegendInterpreter', 'latex')

ax=gca;

if ~isempty(ax.Colorbar)
    ax.Colorbar.Label.Interpreter = 'latex';
    ax.Colorbar.Title.Interpreter = 'latex';
    ax.Colorbar.Label.FontSize = fontS2;
    ax.Colorbar.Title.FontSize = fontS2;
    ax.Colorbar.FontName = 'CMU Serif';
    ax.Colorbar.TickLabelInterpreter='latex';
end

ax.Title.FontSize=fontS2;

set(gca,'FontSize',fontS)
set(gca,'Linewidth',0.5)

set(get(gca,'xlabel'),'Fontsize',fontS2)
set(get(gca,'ylabel'),'Fontsize',fontS2)
set(get(gca,'zlabel'),'Fontsize',fontS2)
set(get(gca,'xlabel'),'Interpreter','latex')
set(get(gca,'ylabel'),'Interpreter','latex')
set(get(gca,'zlabel'),'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
box on

hLegend = findobj(gcf, 'Type', 'Legend');
if ~isempty(hLegend)
    for i=1:numel(hLegend)
        hLegend(i).FontName='CMU Serif';

        hLegend(i).Interpreter='latex';
        hLegend(i).FontSize=fontS;
    end
end

end
