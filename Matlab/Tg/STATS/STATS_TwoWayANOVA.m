function R = STATS_TwoWayANOVA(y, f1, f2, varargin)
% STATS_TwoWayANOVA  Two-way between-subjects ANOVA (Type III) via anovan.
%
% Example call:
%   R = STATS_TwoWayANOVA(d4.Procedural, d4.Genotype, d4.AgeGroup, ...
%                         'Names', {'Genotype','AgeGroup'});
%
% Standard factorial between-subjects ANOVA with the two main effects and
% their interaction. Type III sums of squares (each effect adjusted for all
% others), which is the appropriate choice for the unbalanced cell counts
% here (n differs by Genotype x AgeGroup). H0 for each row: that effect has
% no influence on the mean of y. Use this for single-timepoint designs (e.g.
% Day-4 only) where there is no repeated-measures factor -- do NOT use it
% when a within-subject factor (Day) is present.
%
% INPUT
%   y    numeric response vector (one value per subject)
%   f1   grouping vector for factor 1 (categorical / cellstr / numeric)
%   f2   grouping vector for factor 2
%   varargin: 'Names'({'F1','F2'})  effect labels for output
%
% OUTPUT
%   R.effects  table: Effect, F, DF1, DF2, pValue, partial_eta2
%   R.method   description string
%   R.n        total subjects used (after dropping NaN y)
%   R.tbl      native anovan cell table
%
% partial_eta2 = SS_effect / (SS_effect + SS_error).
% SS 2026

p.Names = {'F1','F2'};
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end

ok = ~isnan(y);
y=y(ok); f1=f1(ok); f2=f2(ok);
[~, tbl] = anovan(y, {f1, f2}, 'model','interaction', 'sstype',3, ...
    'varnames', p.Names, 'display','off');

hdr = tbl(1,:);
cSS = find(contains(hdr,'Sum Sq'),1);
cDF = find(strcmp(hdr,'d.f.'),1);
cF  = find(strcmp(hdr,'F'),1);
cP  = find(contains(hdr,'Prob'),1);
src = string(tbl(2:end,1));

ssErr = tbl{1+find(src=="Error",1), cSS};
dfErr = tbl{1+find(src=="Error",1), cDF};

isEff = ~ismember(src, ["Error","Total"]);
ei = find(isEff);
nE = numel(ei);
Effect = strings(nE,1); F=nan(nE,1); DF1=nan(nE,1); pv=nan(nE,1); etap2=nan(nE,1);
for k = 1:nE                                     % LOOP: 3 effect rows
    r = 1 + ei(k);
    Effect(k) = src(ei(k));
    ss   = tbl{r,cSS};
    F(k) = tbl{r,cF};
    DF1(k)= tbl{r,cDF};
    pv(k)= tbl{r,cP};
    if (ss+ssErr) > 0, etap2(k) = ss/(ss+ssErr); end
end
Effect = replace(Effect, "*", ":");              % 'Genotype*AgeGroup' -> 'Genotype:AgeGroup'

R.effects = table(Effect, F, DF1, repmat(dfErr,nE,1), pv, etap2, ...
    'VariableNames', {'Effect','F','DF1','DF2','pValue','partial_eta2'});
R.method  = 'Two-way between-subjects ANOVA (anovan, Type III SS)';
R.n       = numel(y);
R.tbl     = tbl;
end