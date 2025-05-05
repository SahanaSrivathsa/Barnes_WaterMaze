function saveTablePNG_APA(T, outPNG, varargin)


% ─── options ────────────────────────────────────────────────────────────
p = inputParser;
p.addParameter('Title','',@ischar);
p.addParameter('FontSize',11,@isscalar);           % body font
p.addParameter('Resolution',300,@isscalar);
p.parse(varargin{:});   opt = p.Results;

% ─── 1) round numeric vars & convert to char with two decimals ──────────
for v = 1:width(T)
    if isnumeric(T{1,v})
        T.(v) = arrayfun(@(x) sprintf('%.2f', round(x,4)), T.(v), 'uni',0);
    end
end

% ─── 2) blank duplicate Strategy entries so it appears once per block ───
if ismember('Strategy', T.Properties.VariableNames)
    s = string(T.Strategy);
    T.Strategy([false ; s(2:end)==s(1:end-1)]) = "";
end

% ─── 3) prepare uitable-safe cell array (string → char) ────────────────
raw   = table2cell(T);
isStr = cellfun(@isstring, raw);
raw(isStr) = cellfun(@char, raw(isStr), 'uni',0);

% ─── 4) work out one uniform column width (in pixels) ───────────────────
%       crude but reliable: (#characters in widest entry) * 7 px + padding
wc = cellfun(@numel, [T.Properties.VariableNames ; raw(:)]);
colWidth = num2cell( (max(wc)+2)*7 * ones(1,width(T)) );   % same for all

% ─── 5) figure & uitable ────────────────────────────────────────────────
rowH  = 24;                                % px per row
figH  = rowH*(height(T)+2);                % + title
figW  = (max(wc)+2)*7*width(T)+40;         % crude total width
fig   = figure('Visible','on','Units','pixels',...
               'Position',[100 100 figW figH],'Color','w');

hdr = strcat("<html><b>", T.Properties.VariableNames, "</b></html>");

uitable(fig,'Data',raw,...
              'ColumnName',cellstr(hdr),...
              'RowName','',...
              'Units','normalized',...
              'Position',[0 0 1 0.92],...
              'FontSize',opt.FontSize,...
              'ColumnWidth',colWidth);

if ~isempty(opt.Title)
    annotation(fig,'textbox',[0 0.92 1 0.08],...
               'String',opt.Title,...
               'HorizontalAlignment','center',...
               'FontWeight','bold',...
               'EdgeColor','none',...
               'FontSize',14);          % ← title size 14
end

drawnow;
im = getframe(fig);
imwrite(im.cdata, outPNG);
close(fig);
end
