function pubify_figure_axis_robust(fontSizeX, fontSizeY)
% Makes figures publication-ready
    ax = gca;
    ax.FontSize = max(fontSizeX, fontSizeY);
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLabel.FontSize = fontSizeX;
    ax.YLabel.FontSize = fontSizeY;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
    if ~isempty(ax.Title)
        ax.Title.FontSize = max(fontSizeX, fontSizeY) + 2;
        ax.Title.FontWeight = 'bold';
    end
end