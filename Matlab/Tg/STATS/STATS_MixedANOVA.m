function R = STATS_MixedANOVA(L, valueVar, withinVars, betweenVars, subjectVar)
% STATS_MixedANOVA  Mixed-design repeated-measures ANOVA via fitrm/ranova.
%

% SS 2026

%% Validate
withinVars  = cellstr(withinVars);
betweenVars = cellstr(betweenVars);
need = [{valueVar, subjectVar}, withinVars, betweenVars];
have = L.Properties.VariableNames;
miss = need(~ismember(need, have));
if ~isempty(miss)
    error('STATS_MixedANOVA:missingCols','L missing column(s): %s', strjoin(miss,', '));
end
if isempty(withinVars)
    error('STATS_MixedANOVA:noWithin','At least one within-subject factor is required.');
end
nw = numel(withinVars);

%% Within-factor levels (categorical order if categorical, else sorted unique)
wlev = cell(1,nw);
for j = 1:nw                                     % LOOP: <=2 within factors
    col = L.(withinVars{j});
    if iscategorical(col)
        wlev{j} = categories(removecats(col));
    else
        wlev{j} = num2cell(unique(col));         % numeric/other -> ascending
    end
end

%% Full crossing of within levels -> combo table (defines wide-column order)
% Day varies fastest within each level of the next factor.
sz = cellfun(@numel, wlev);
nCombo = prod(sz);
combo = table();
for j = 1:nw
    rep = repmat(reshape(repmat((1:sz(j)), prod(sz(1:j-1)), 1), [], 1), prod(sz(j+1:end)), 1);
    lv  = wlev{j};
    if isnumeric(lv{1}) || islogical(lv{1})
        combo.(withinVars{j}) = cell2mat(lv(rep));
    else
        combo.(withinVars{j}) = categorical(lv(rep), lv);
    end
end

%% Build wide response matrix Y (subjects x combos), vectorised match by key
subs = unique(L.(subjectVar), 'stable');
ns   = numel(subs);
keyL = string(L.(subjectVar));
for j = 1:nw, keyL = keyL + "|" + string(L.(withinVars{j})); end
Y = nan(ns, nCombo);
colNames = "y_" + (1:nCombo);
for c = 1:nCombo                                 % LOOP: <=12 combos
    ck = string(subs);
    for j = 1:nw, ck = ck + "|" + string(combo.(withinVars{j})(c)); end
    [tf, loc] = ismember(ck, keyL);
    Y(tf, c)  = L.(valueVar)(loc(tf));
end

%% Between-subject values per subject (assumed constant within subject)
[~, firstIdx] = unique(L.(subjectVar), 'stable');
W = array2table(Y, 'VariableNames', cellstr(colNames));
for b = 1:numel(betweenVars)
    W.(betweenVars{b}) = categorical(L.(betweenVars{b})(firstIdx));
end

%% Drop incomplete subjects (fitrm needs complete within structure)
bad = any(isnan(Y), 2);
if any(bad)
    fprintf('[WARN] STATS_MixedANOVA: dropping %d/%d subjects with incomplete within-cells\n', nnz(bad), ns);
    W = W(~bad, :);
end

%% Fit fitrm + ranova
respRange   = sprintf('%s-%s', char(colNames(1)), char(colNames(end)));
betweenModel = '1';
if ~isempty(betweenVars), betweenModel = strjoin(betweenVars, '*'); end
rm = fitrm(W, sprintf('%s ~ %s', respRange, betweenModel), 'WithinDesign', combo);

withinModel = strjoin(withinVars, '*');
rv  = ranova(rm, 'WithinModel', withinModel);
bwn = anova(rm);                                 % between-subjects tests

%% Tidy effects table with partial eta^2 = SS / (SS + SS_error)
rn    = string(rv.Properties.RowNames);
isErr = startsWith(rn, "Error");
ssAll = rv.SumSq;
betweenErrSS = ssAll(rn == "Error");

keep = ~isErr & rn ~= "(Intercept)";
ke   = find(keep);
nE   = numel(ke);
etap2 = nan(nE,1); df2 = nan(nE,1);
for k = 1:nE                                     % LOOP: small effects table
    toks   = split(rn(ke(k)), ":");
    inWith = toks(ismember(toks, string(withinVars)));   % within factors in this effect
    if isempty(inWith)
        errSS  = betweenErrSS;                   % pure between effect
        df2(k) = rv.DF(rn == "Error");
    else
        ord    = arrayfun(@(x) find(string(withinVars)==x,1), inWith);
        [~,o]  = sort(ord);
        errNm  = "Error(" + strjoin(inWith(o), ":") + ")";
        errSS  = ssAll(rn == errNm);
        df2(k) = rv.DF(rn == errNm);
    end
    if ~isempty(errSS) && (ssAll(ke(k)) + errSS) > 0
        etap2(k) = ssAll(ke(k)) / (ssAll(ke(k)) + errSS);
    end
end

effName = regexprep(rn(ke), '^\(Intercept\):', '');   % '(Intercept):Day' -> 'Day'
R.effects = table(effName, rv.F(ke), rv.DF(ke), df2, rv.pValue(ke), etap2, ...
    'VariableNames', {'Effect','F','DF1','DF2','pValue','partial_eta2'});

%% n per between-cell
if ~isempty(betweenVars)
    R.n = groupsummary(W, betweenVars);
else
    R.n = table(height(W), 'VariableNames', {'GroupCount'});
end
R.n_total     = height(W);
R.method      = 'Mixed-design RM-ANOVA (fitrm/ranova, Type III SS)';
R.ranova_raw  = rv;
R.between_raw = bwn;
R.rm          = rm;
R.within      = withinVars;
R.between     = betweenVars;
end