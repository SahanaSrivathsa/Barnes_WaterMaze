function R = STATS_MWM_LME(T, varName, varargin)
% R = STATS_MWM_LME(T, 'Platform_CIPL')  % Genotype*Sex (default)
% R = STATS_MWM_LME(T, 'thigmotaxis', 'Factors', "Genotype")  % genotype-only (single-sex data)
%
% Fit a linear mixed-effects model (repeated days per rat fitlme 
%  or a general linear model for Morris watermaze variable
%
% Between-subject fixed effects are whatever is listed in 'Factors'
% Default-Genotype*Sex (+ Day) for repeated data
% 
% ONLY MALES:P
% Input 'Factors', "Genotype" for a single-sex cohort.
%

% Repeated daily measures per rat are correlated
% random intercept gets between-rat variance so the fixed
% effects are tested against the correct (within-rat) error. Gaussian errors
% are assumed; residual normality is checked .
%
% SS 2026

% ---- defaults ----
p.Factors         = ["Genotype","Sex"];
p.GenoRef         = 'WT';
p.SexRef          = 'F';
p.DropSingleLevel = true;
for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end
facReq = string(p.Factors);

% ---- validate required columns ----
req  = ["Animal", string(varName), facReq];
have = string(T.Properties.VariableNames);
miss = req(~ismember(req, have));
if ~isempty(miss)
    error('STATS_MWM_LME:missingCols', ...
        'T missing column(s): %s', strjoin(cellstr(miss), ', '));
end

% ---- drop NaN DV ----
T = T(~isnan(T.(varName)), :);

% ---- decide which factors to keep (>= 2 levels), drop missing rows for them ----
refMap = containers.Map({'Genotype','Sex'}, {string(p.GenoRef), string(p.SexRef)});
facUse = strings(0,1);
for f = facReq(:)'                                  % LOOP: tiny factor list
    levs = unique(string(T.(f)));
    levs = levs(~ismissing(levs) & levs ~= "");
    if numel(levs) < 2
        if p.DropSingleLevel
            warning('STATS_MWM_LME:singleLevel', ...
                'Factor "%s" has <2 levels for %s; dropped from model.', f, varName);
            continue
        else
            error('STATS_MWM_LME:singleLevel', ...
                'Factor "%s" has <2 levels for %s.', f, varName);
        end
    end
    facUse(end+1,1) = f;                              %#ok<AGROW> tiny list
end

% drop rows with missing levels of the retained factors
for f = facUse(:)'
    T = T(~ismissing(T.(f)), :);
end
if height(T) < 4
    error('STATS_MWM_LME:tooFewRows', ...
        'Only %d valid rows for %s after dropping NaN/missing.', height(T), varName);
end
sd_dv = std(T.(varName), 'omitnan');

% ---- set categorical ordering (reference level first) for retained factors ----
for f = facUse(:)'
    levs = unique(string(T.(f)));
    if isKey(refMap, char(f)), ref = refMap(char(f)); else, ref = levs(1); end
    if ~ismember(ref, levs)
        warning('STATS_MWM_LME:refMissing', ...
            'Reference %s not found for %s; using %s.', ref, f, levs(1));
        ref = levs(1);
    end
    order   = [ref; levs(levs ~= ref)];
    T.(f)   = categorical(string(T.(f)), cellstr(order));
end

% ---- build the between-subject part of the formula ----
if isempty(facUse)
    between = '1';                                   % intercept-only
else
    between = strjoin(cellstr(facUse), '*');         % e.g. 'Genotype*Sex' or 'Genotype'
end

% ---- repeated (Day) vs single-trial ----
nDayPerRat = max(accumarray(findgroups(T.Animal), ones(height(T),1)));
isRepeated = (nDayPerRat > 1) && ismember('Day', T.Properties.VariableNames);

% ---- FIT MODEL ----
if isRepeated
    T.Day_c = T.Day - mean(T.Day,'omitnan');
    rhs = [between ' + Day_c'];
    if ~isempty(facUse)                              % add each factor x Day_c (2-way)
        rhs = [rhs ' + ' strjoin(strcat(cellstr(facUse), ':Day_c'), ' + ')];
    end
    rhs       = [rhs ' + (1|Animal)'];
    formula   = sprintf('%s ~ %s', varName, rhs);
    mdl       = fitlme(T, formula);
    termNames = string(mdl.Coefficients.Name);
    beta      = mdl.Coefficients.Estimate;
    se        = mdl.Coefficients.SE;
    tstat     = mdl.Coefficients.tStat;
    df        = mdl.Coefficients.DF;
    pval      = mdl.Coefficients.pValue;
    ci_mat    = [mdl.Coefficients.Lower, mdl.Coefficients.Upper];
    resid     = mdl.residuals;
    modelType = 'LME_Repeated';
    method    = ['fitlme: ' formula];
else
    formula   = sprintf('%s ~ %s', varName, between);
    mdl       = fitlm(T, formula);
    termNames = string(mdl.CoefficientNames)';
    beta      = mdl.Coefficients.Estimate;
    se        = mdl.Coefficients.SE;
    tstat     = mdl.Coefficients.tStat;
    df        = repmat(mdl.DFE, numel(termNames), 1);
    pval      = mdl.Coefficients.pValue;
    ci_mat    = coefCI(mdl);
    resid     = mdl.Residuals.Raw;
    modelType = 'LM_SingleTrial';
    method    = ['fitlm: ' formula];
end

% ---- effect size (standardised beta) ----
stdBeta = beta ./ sd_dv;

fixedFX = table(termNames, beta, se, tstat, df, pval, ci_mat(:,1), ci_mat(:,2), stdBeta, ...
    'VariableNames', {'Term','Estimate','SE','tStat','DF','pValue','Lower','Upper','std_beta'});

% ---- residual normality ----
if numel(resid) >= 8
    [~, norm_p] = lillietest(resid);
else
    norm_p = NaN;
end
norm_flag = 'OK';
if ~isnan(norm_p) && norm_p < 0.05
    norm_flag = 'NON-NORMAL (consider log transform)';
    warning('STATS_MWM_LME:residNorm', ...
        '%s residuals non-normal (Lilliefors p=%.3f) for %s.', modelType, norm_p, varName);
end

% ---- output ----
R.fixedEffects = fixedFX;
R.n_obs        = height(T);
R.n_rats       = numel(unique(T.Animal));
R.residNorm_p  = norm_p;
R.norm_flag    = norm_flag;
R.model_type   = modelType;
R.formula      = formula;
R.method       = method;
R.factors_used = facUse;
R.n            = R.n_obs;
R.p            = pval;
R.stat         = tstat;
R.ci           = ci_mat;
R.effect_size  = stdBeta;
end