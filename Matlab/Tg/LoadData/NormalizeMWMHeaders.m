function [T, renamed] = NormalizeMWMHeaders(T, varargin)
%
% Rename the columns of a Morris-Water-Maze sheet to canonical names, so that
% files exported with different header spellings concatenate cleanly. C
alias = containers.Map( ...
  {'test','animal','trial','stage','apparatus','duration','distance','meanspeed', ...
   'pathefficiency','platformcipl','platformentries', ...
   'q1time','q2time','q3time','q4time','platformpathefficiencytoentry'}, ...
  {'Test','Animal','Trial','Stage','Apparatus','Duration','Distance','MeanSpeed', ...
   'PathEfficiency','Platform_CIPL','Platform_Entries', ...
   'Q1_Time','Q2_Time','Q3_Time','Q4_Time','Platform_PathEfficiencyToEntry'});
for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'Alias'), alias = varargin{k+1}; end
end

raw   = string(T.Properties.VariableNames);
nCol  = numel(raw);
keys  = strings(nCol,1);
canon = raw(:);                          % default: unchanged
matched = false(nCol,1);

for c = 1:nCol                           % LOOP: per-column header lookup (small, not vectorizable)
    key     = regexprep(lower(raw(c)), '\s*\([^)]*\)', '');   % strip unit parentheticals
    key     = regexprep(key, '[^a-z0-9]', '');                % drop punctuation/spaces
    keys(c) = key;
    if isKey(alias, char(key))
        canon(c)   = string(alias(char(key)));
        matched(c) = true;
    end
end

T.Properties.VariableNames = cellstr(canon);
renamed = table(raw(:), keys, canon, matched, ...
    'VariableNames', {'Raw','Key','Canonical','Matched'});
end