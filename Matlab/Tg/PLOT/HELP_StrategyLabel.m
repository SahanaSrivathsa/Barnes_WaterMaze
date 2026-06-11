function lab = HELP_StrategyLabel(name)
% lab = HELP_StrategyLabel(name)
%
% Map a strategy/variable COLUMN name to its display label. Data columns keep
% their original names (circling, chaining); only the human-facing label
% changes. Unmapped names pass through unchanged.
%   circling -> Spiralling
%   chaining -> Ring search
%
% SS 2026

name = char(string(name));
switch lower(name)
    case 'circling', lab = 'Spiralling';
    case 'chaining', lab = 'Ring search';
    otherwise,       lab = name;
end
end