function stylePaperFigure(fig)
%STYLEPAPERFIGURE Apply the paper figure style used by result scripts.
if nargin < 1 || isempty(fig)
    fig = gcf;
end

axesList = findall(fig, 'Type', 'axes');
for iax = 1:numel(axesList)
    ax = axesList(iax);
    set(ax, ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'LineWidth', 1.2, ...
        'Box', 'on');

    styleLabel(ax.XLabel);
    styleLabel(ax.YLabel);
    styleLabel(ax.ZLabel);
    styleLabel(ax.Title);

    try
        ax.XAxis.FontSize = 14;
        ax.XAxis.FontWeight = 'bold';
        ax.YAxis.FontSize = 14;
        ax.YAxis.FontWeight = 'bold';
    catch
        % Older MATLAB releases do not expose axis rulers consistently.
    end
end

legendList = findall(fig, 'Type', 'legend');
for ilg = 1:numel(legendList)
    lgd = legendList(ilg);
    entryCount = numLegendEntries(lgd);
    set(lgd, ...
        'FontSize', 13, ...
        'FontWeight', 'bold', ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off');
    try
        if entryCount <= 3
            lgd.NumColumns = max(1, entryCount);
        elseif entryCount <= 6
            lgd.NumColumns = 3;
        else
            lgd.NumColumns = 4;
        end
    catch
    end
end
end

function styleLabel(h)
if isempty(h) || ~isvalid(h)
    return;
end
set(h, 'FontSize', 16, 'FontWeight', 'bold');
end

function n = numLegendEntries(lgd)
labels = lgd.String;
if ischar(labels)
    n = size(labels, 1);
else
    n = numel(labels);
end
end
