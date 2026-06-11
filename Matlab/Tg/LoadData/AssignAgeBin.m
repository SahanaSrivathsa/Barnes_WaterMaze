function ageGroup = AssignAgeBin(age_mo, edges, labels)
%
% Assign a continuous age (in MONTHS) to a named month-bin. Binning is done
% with DISCRETIZE, so bins are half-open [edges(k), edges(k+1)). P

labels = string(labels(:))';
validateattributes(edges, {'numeric'}, {'vector', 'increasing'}, mfilename, 'edges');
if numel(labels) ~= numel(edges) - 1
    error('HELP_AssignAgeBin:sizeMismatch', ...
        'labels must have numel(edges)-1 (=%d) elements, got %d.', numel(edges)-1, numel(labels));
end

binIdx  = discretize(age_mo(:), edges);          % NaN where out of range / NaN age
inRange = ~isnan(binIdx);

lab            = repmat("", numel(age_mo), 1);
lab(inRange)   = labels(binIdx(inRange));
ageGroup       = categorical(lab, labels, 'Ordinal', true); 

if any(~inRange)
    warning('HELP_AssignAgeBin:outOfRange', ...
        '%d age value(s) fell outside [%.1f, %.1f) (or were NaN) and are <undefined>.', ...
        nnz(~inRange), edges(1), edges(end));
end
end