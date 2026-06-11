function S = PLOT_Tg_DeltaByAge(T, varName, varargin)
% PLOT_Tg_DeltaByAge  Day4-Day3 change in strategy probability, WT vs APP by age.
%
% Example call:
%   S = PLOT_Tg_DeltaByAge(T, 'Procedural', 'AgeCol','AgeGroup', ...
%           'Colors', genoCols, 'SaveDir', dir_out);
%
% Per-rat change score Delta = P(Day d2) - P(Day d1) for varName, summarised
% as a grouped bar (mean +/- SEM) of WT vs APP within each age bin, with
% per-rat dots and a y=0 reference line. Delta makes directionality explicit:
% Delta<0 = strategy dropped from d1->d2, Delta>0 = rose.
%
% Two test families, each Bonferroni-corrected within family:
%   (A) DIRECTIONALITY  one-sample t of Delta vs 0, per Genotype x age
%       (identical to the paired Day d1-vs-d2 test). H0: mean change = 0.
%       Corrected over all valid cells (default m = #cells). Significant
%       cells get a star drawn at the bar.
%   (B) INTERACTION     WT vs APP on Delta, per age (Welch t2 by default).
%       H0: the d1->d2 change is the same in both genotypes. This equals the
%       Genotype x Day interaction of a 2-day mixed ANOVA; opposite-sign means
%       with a significant test = crossover. Corrected over ages (default m =
%       #ages). Drawn via sigstar using only the CORRECTED p-values.
%
% Welch is used for (B) because WT/APP variances and n differ; flip to
% 'ranksum' for a rank test on small/old cells.
%
% INPUT
%   T        per-rat-per-day strategy table (Animal, Genotype, AgeCol, Day, varName)
%   varName  strategy-group column to analyse (e.g. 'Procedural')
%   varargin 'Days'([3 4]), 'AgeCol'('AgeGroup'), 'AgeLevels'(["Young","Mid","Old"]),
%            'Genotypes'(["WT","APP"]), 'Colors'(geno), 'Test'('ttest2'|'ranksum'),
%            'Correction'('bonferroni'|'none'), 'ShowDirStars'(true),
%            'YLabel', 'SaveDir'.
%
% OUTPUT
%   S        long stats table: comparison ('Delta_vs0' | 'WT_vs_APP'), means,
%            stat, df, p_raw, p_corr, m, test, note (direction / opposite).
%
% USES: HELP_DrawGroupedBars, sigstar, HELP_StrategyLabel.
% SS 2026

p.Days=[3 4]; p.AgeCol='AgeGroup'; p.AgeLevels=["Young","Mid","Old"];
p.Genotypes=["WT","APP"]; p.Colors=[0.50 0.70 1.00; 0.00 0.20 0.70];
p.Test='ttest2'; p.Correction='bonferroni'; p.ShowDirStars=true;
p.YLabel=''; p.SaveDir='';
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end

req={'Animal','Day','Genotype',p.AgeCol,char(varName)};
have=T.Properties.VariableNames; miss=req(~ismember(req,have));
if ~isempty(miss)
    error('PLOT_Tg_DeltaByAge:missingCols','T missing column(s): %s',strjoin(miss,', '));
end
d1=p.Days(1); d2=p.Days(2);
lab = char(varName); if exist('HELP_StrategyLabel','file'), lab=HELP_StrategyLabel(varName); end
if isempty(p.YLabel), p.YLabel=sprintf('\\DeltaP(%s)  (Day %d - Day %d)', lab, d2, d1); end

%% Per-rat change score Delta = P(d2) - P(d1)
Td = T(ismember(T.Day,[d1 d2]), {'Animal','Genotype',p.AgeCol,'Day',char(varName)});
Td.DayLab = categorical("D"+string(Td.Day));
U = unstack(Td(:,{'Animal','Genotype',p.AgeCol,'DayLab',char(varName)}), char(varName), 'DayLab');
colA = char("D"+d1); colB = char("D"+d2);
U.Delta = U.(colB) - U.(colA);
U = U(~isnan(U.Delta), :);

G = string(p.Genotypes); nG=numel(G);
ages = string(p.AgeLevels);
ages = ages(ismember(ages, unique(string(U.(p.AgeCol)))));   % keep populated, in order
ages = ages(:);                                              % column (table-safe indexing)
nA = numel(ages);
if nA==0, warning('PLOT_Tg_DeltaByAge:noData','No populated age levels.'); S=table(); return; end

%% Collect per-cell Delta vectors
vals = cell(nA,nG);
for a=1:nA                                       % LOOP: ages x genotypes (small)
    for g=1:nG
        vals{a,g} = U.Delta(U.(p.AgeCol)==ages(a) & string(U.Genotype)==G(g));
    end
end

%% (A) directionality: one-sample Delta vs 0 per cell
tA=nan(nA,nG); dfA=nan(nA,nG); pA=nan(nA,nG); mA=nan(nA,nG); nA_=zeros(nA,nG);
for a=1:nA
    for g=1:nG
        v=vals{a,g}; nA_(a,g)=numel(v);
        if numel(v)>1
            mA(a,g)=mean(v);
            [~,pp,~,st]=ttest(v); tA(a,g)=st.tstat; dfA(a,g)=st.df; pA(a,g)=pp;
        elseif numel(v)==1, mA(a,g)=v;
        end
    end
end

%% (B) interaction: WT vs APP on Delta per age
tB=nan(nA,1); dfB=nan(nA,1); pB=nan(nA,1); testB=strings(nA,1);
for a=1:nA
    w=vals{a,1}; c=vals{a,2};
    if numel(w)>1 && numel(c)>1
        if strcmpi(p.Test,'ranksum')
            [pp,~,st]=ranksum(w,c); testB(a)="ranksum"; if isfield(st,'zval'), tB(a)=st.zval; end
        else
            [~,pp,~,st]=ttest2(w,c,'Vartype','unequal'); tB(a)=st.tstat; dfB(a)=st.df; testB(a)="ttest2_Welch";
        end
        pB(a)=pp;
    end
end

%% Bonferroni within each family
switch lower(p.Correction)
    case 'bonferroni'
        mA_n=nnz(~isnan(pA)); pAc=min(pA*mA_n,1);
        mB_n=nnz(~isnan(pB)); pBc=min(pB*mB_n,1);
    case 'none'
        mA_n=NaN; pAc=pA; mB_n=NaN; pBc=pB;
    otherwise, error('PLOT_Tg_DeltaByAge:badCorr','Unknown Correction: %s',p.Correction);
end

%% Draw
f=figure('Position',[80 80 520+120*nA 460]); ax=gca;
[hBar,ctr]=HELP_DrawGroupedBars(ax,vals,p.Colors,G,cellstr(ages),0.85,0.12);
yline(ax,0,'k-','LineWidth',0.75,'HandleVisibility','off');
ylabel(ax,p.YLabel); xlabel(ax,'Age bin');
title(ax,sprintf('%s  -  change Day %d -  (WT vs APP)',lab,d1,d2),'FontWeight','bold');
legend(hBar,cellstr(G),'Location','northeastoutside');
if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(13,13); end

% (B) between-genotype sigstar with corrected p
sigPairs={}; sigP=[];
for a=1:nA
    if ~isnan(pBc(a)) && pBc(a)<=0.05
        sigPairs{end+1}=[ctr(a,1),ctr(a,2)]; sigP(end+1)=pBc(a); %#ok<AGROW>
    end
end
if ~isempty(sigPairs) && exist('sigstar','file'), sigstar(sigPairs,sigP); end

% (A) directionality stars at each bar (corrected)
if p.ShowDirStars
    yr=diff(ylim); 
    for a=1:nA
        for g=1:nG
            if ~isnan(pAc(a,g)) && pAc(a,g)<=0.05
                if pAc(a,g)<=1e-3, s='***'; elseif pAc(a,g)<=1e-2, s='**'; else, s='*'; end
                yoff = sign(mA(a,g))*0.03*yr; if mA(a,g)==0, yoff=0.03*yr; end
                va = 'bottom'; if mA(a,g)<0, va='top'; end
                text(ax,ctr(a,g), mA(a,g)+yoff, s, 'HorizontalAlignment','center', ...
                    'VerticalAlignment',va, 'FontSize',14, 'Color',p.Colors(g,:)*0.6);
            end
        end
    end
end

%% Stats table (long): directionality rows + interaction rows
Sa = table();
[gg,aa]=ndgrid(1:nG,1:nA);
Sa.Variable   = repmat(string(lab),nG*nA,1);
Sa.AgeGroup   = ages(aa(:));
Sa.comparison = repmat("Delta_vs0",nG*nA,1);
Sa.Genotype   = G(gg(:))';
Sa.n1         = nA_(sub2ind([nA nG],aa(:),gg(:)));
Sa.n2         = nan(nG*nA,1);
Sa.mean1      = mA(sub2ind([nA nG],aa(:),gg(:)));
Sa.mean2      = zeros(nG*nA,1);
Sa.stat       = tA(sub2ind([nA nG],aa(:),gg(:)));
Sa.df         = dfA(sub2ind([nA nG],aa(:),gg(:)));
Sa.p_raw      = pA(sub2ind([nA nG],aa(:),gg(:)));
Sa.p_corr     = pAc(sub2ind([nA nG],aa(:),gg(:)));
Sa.m          = repmat(mA_n,nG*nA,1);
Sa.test       = repmat("ttest_vs0",nG*nA,1);
Sa.note       = repmat("",nG*nA,1);
Sa.note(Sa.mean1<0) = "down"; Sa.note(Sa.mean1>0) = "up";

Sb = table();
Sb.Variable   = repmat(string(lab),nA,1);
Sb.AgeGroup   = ages(:);
Sb.comparison = repmat("WT_vs_APP",nA,1);
Sb.Genotype   = repmat("WT-APP",nA,1);
Sb.n1         = nA_(:,1); Sb.n2 = nA_(:,2);
Sb.mean1      = mA(:,1);  Sb.mean2 = mA(:,2);
Sb.stat       = tB; Sb.df = dfB; Sb.p_raw = pB; Sb.p_corr = pBc;
Sb.m          = repmat(mB_n,nA,1); Sb.test = testB;
Sb.note       = repmat("",nA,1);
Sb.note(mA(:,1).*mA(:,2)<0) = "opposite";
S = [Sa; Sb];
S.correction = repmat(string(p.Correction),height(S),1);

if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    tok = strrep(lab,' ','_');
    saveas(f,fullfile(p.SaveDir,sprintf('Delta_Day%dv%d_%s_%s.png',d1,d2,tok,p.AgeCol)));
    close(f);
end
end