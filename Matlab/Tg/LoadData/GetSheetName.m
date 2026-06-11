function sheetName = GetSheetName(fname, aliases)
% Gets all sheet names in excel spatial file

aliases = cellstr(aliases);
 
[~, rawSheets] = xlsfinfo(fname);
if isempty(rawSheets)
    error('GetSheetName:unreadable', 'Could not read sheet list from %s', fname);
end
% Normalize to a cell of chars regardless of MATLAB version's return type.
sheets    = cellstr(rawSheets);
sheetsLow = lower(string(sheets));
 
for k = 1:numel(aliases)             % LOOP: small fixed candidate list
    if any(sheetsLow == lower(string(aliases{k})))
        sheetName = aliases{k};
        return
    end
end
 
error('GetSheetName:notFound', ...
    'No matching sheet in %s.\n  Tried   : {%s}\n  Present : {%s}', ...
    fname, strjoin(aliases, ', '), strjoin(sheets, ', '));
end