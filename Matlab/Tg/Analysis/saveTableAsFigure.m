function saveTableAsFigure(t, fname, varargin)
% saveTablePNG   Snapshot a table into a PNG (UI-safe, row names kept).
%
%   saveTablePNG(t, 'anovaYoung.png', 'Title','ANOVA – Young');
%
% Optional name-value:
%   'Title'       text shown at top of figure
%   'FontSize'    uitable font size  (default 12)
%   'RowName'     'numbered' | [] | cellstr (default t.RowNames or numbered)

% ---- parse inputs ----
p = inputParser;
p.addParameter('Title', '', @ischar);
p.addParameter('FontSize', 12, @isscalar);
p.addParameter('RowName',  'auto');
p.parse(varargin{:});
opt = p.Results;

if strcmpi(opt.RowName,'auto')
    if ~isempty(t.Properties.RowNames)
        rowNameOpt = t.Properties.RowNames;
    else
        rowNameOpt = 'numbered';
    end
else
    rowNameOpt = opt.RowName;
end

% ---- build invisible figure with uitable ----
rowH   = 22;                           % pixels per row
figH   = rowH*(height(t)+2);           % +2 for padding/title
figW   = 900;
fig    = figure('Visible','on', ...
                'Units','pixels', ...
                'Position',[100 100 figW figH], ...
                'Color','w');

% uitable occupying full figure
uitable(fig, ...
    'Data'      , table2cell(t), ...
    'ColumnName', t.Properties.VariableNames, ...
    'RowName'   , rowNameOpt, ...
    'Units'     , 'normalized', ...
    'Position'  , [0 0 1 0.95], ...
    'FontSize'  , opt.FontSize);

if ~isempty(opt.Title)
    annotation(fig,'textbox',[0 0.95 1 0.05], ...
               'String', opt.Title, ...
               'HorizontalAlignment','center', ...
               'FontWeight','bold', ...
               'EdgeColor','none', ...
               'FontSize',opt.FontSize+2);
end

% ---- make sure UI layer is painted, then snapshot ----
drawnow;                       % critical!
im = getframe(fig);            % includes uitable w/o warnings
imwrite(im.cdata, fname);      % e.g. 'anovaYoung.png'
close(fig);
end
