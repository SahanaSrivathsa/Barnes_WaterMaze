function S = PLOT_Tg_Day4Groups_Panel(T, binCol, varargin)
% PLOT_Tg_Day4Groups_Panel  Day-4 WT vs APP per strategy group, corrected stars.
%
% Example call:
%   S = PLOT_Tg_Day4Groups_Panel(T, 'AgeGroup', 'Day', 4, ...
%           'Colors', genoCols, 'Correction', 'bonferroni', 'SaveDir', dir_out);
%
% One slide-ready 1 x 3 figure (Platform-Independent / Procedural /
% Allocentric). Each panel = WT vs APP grouped bar of Day-`Day` group
% probability across the age bins of binCol, with per-cell WT-vs-APP test.
% All cell p-values are collected across the three panels, corrected
% together (Bonferroni by default, m = number of valid tests), and only the
% CORRECTED p-values are passed to sigstar so stars reflect the corrected
% significance. Raw and corrected p are both returned in S.
%
% varargin: 'Day'(4), 'Groups', 'GroupLabels', 'Genotypes'(["WT","APP"]),
%           'Colors'(geno), 'Test'('ttest2'|'ranksum'),
%           'Correction'('bonferroni'|'none'), 'CorrectN'([]=auto),
%           'YLabel', 'SaveDir'.
% SS 2026

p.Day=4;
p.Groups={'PlatformIndependent','Procedural','Allocentric'};
p.GroupLabels={'Platform-Independent','Procedural','Allocentric'};
p.Genotypes=["WT","APP"];
p.Colors=[0.50 0.70 1.00; 0.00 0.20 0.70];
p.Test='ttest2'; p.Correction='bonferroni'; p.CorrectN=[];
p.YLabel='P(strategy group)'; p.SaveDir='';
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end

req=[{'Animal','Day','Genotype',char(binCol)}, p.Groups];
have=T.Properties.VariableNames; miss=req(~ismember(req,have));
if ~isempty(miss)
    error('PLOT_Tg_Day4Groups_Panel:missingCols','T missing column(s): %s',strjoin(miss,', '));
end

Td=T(T.Day==p.Day,:);
bins=categories(removecats(Td.(binCol)));
G=string(p.Genotypes); nB=numel(bins); nG=numel(G); nGrp=numel(p.Groups);

% ---- Pass 1: per-cell values + raw stats (no drawing yet) -------------
valsAll = cell(nGrp,1);            % valsAll{gi}{b,g}
ctrAll  = cell(nGrp,1);            % bar-centre pairs per group
praw  = nan(nGrp,nB); stat=nan(nGrp,nB); df=nan(nGrp,nB); test=strings(nGrp,nB);
mW=nan(nGrp,nB); mA=nan(nGrp,nB); nW=zeros(nGrp,nB); nA=zeros(nGrp,nB);
for gi=1:nGrp                                    % LOOP: 3 strategy groups
    gv=p.Groups{gi};
    rat=groupsummary(Td,{char(binCol),'Genotype','Animal'},'mean',gv);
    ratCol=['mean_' gv];
    vals=cell(nB,nG);
    for b=1:nB
        for g=1:nG
            v=rat.(ratCol)(rat.(binCol)==bins{b} & string(rat.Genotype)==G(g));
            vals{b,g}=v(~isnan(v));
        end
        a=vals{b,1}; c=vals{b,2};
        nW(gi,b)=numel(a); nA(gi,b)=numel(c);
        if ~isempty(a), mW(gi,b)=mean(a); end
        if ~isempty(c), mA(gi,b)=mean(c); end
        if numel(a)>1 && numel(c)>1
            if strcmpi(p.Test,'ranksum')
                [pp,~,st]=ranksum(a,c); test(gi,b)="ranksum";
                if isfield(st,'zval'), stat(gi,b)=st.zval; end
            else
                [~,pp,~,st]=ttest2(a,c,'Vartype','unequal');
                stat(gi,b)=st.tstat; df(gi,b)=st.df; test(gi,b)="ttest2_Welch";
            end
            praw(gi,b)=pp;
        end
    end
    valsAll{gi}=vals;
end

% ---- Multiple-comparison correction across ALL cells together --------
m = p.CorrectN;
if isempty(m), m = nnz(~isnan(praw)); end        % default: # valid tests (= 9 here)
switch lower(p.Correction)
    case 'bonferroni', pcorr = min(praw*m, 1);
    case 'none',       pcorr = praw;
    otherwise, error('PLOT_Tg_Day4Groups_Panel:badCorr','Unknown Correction: %s',p.Correction);
end

% ---- Pass 2: draw bars, then corrected sigstar -----------------------
f=figure('Position',[60 80 1500 460]);
tl=tiledlayout(1,nGrp,'TileSpacing','tight','Padding','tight');
hKeep=gobjects(nG,1);
for gi=1:nGrp                                    % LOOP: one tile per group
    ax=nexttile;
    [hBar,ctr]=HELP_DrawGroupedBars(ax,valsAll{gi},p.Colors,G,bins,0.75,0.12);
    ctrAll{gi}=ctr;
    if gi==1, hKeep=hBar; end
    title(ax,p.GroupLabels{gi},'FontWeight','bold');
    if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(12,12); end

    sigPairs={}; sigP=[];
    for b=1:nB
        if ~isnan(pcorr(gi,b)) && pcorr(gi,b)<=0.05   % only corrected-significant
            sigPairs{end+1}=[ctr(b,1),ctr(b,2)]; sigP(end+1)=pcorr(gi,b); %#ok<AGROW>
        end
    end
    if ~isempty(sigPairs) && exist('sigstar','file'), sigstar(sigPairs,sigP); end
end

xlabel(tl,'Age bin','FontSize',14,'FontWeight','bold');
ylabel(tl,p.YLabel,'FontSize',14,'FontWeight','bold');
title(tl,sprintf('Strategy groups  -  Day %d  (WT vs APP)',p.Day),'FontSize',16,'FontWeight','bold');
legend(hKeep,cellstr(G),'Orientation','horizontal','Location','southoutside');

% ---- Tidy stats table (raw + corrected) ------------------------------
[gg,bb]=ndgrid(1:nGrp,1:nB);
Variable = string(p.GroupLabels(gg(:)))';
AgeBin   = string(bins(bb(:)));
S=table(Variable, repmat(string(binCol),nGrp*nB,1), AgeBin, repmat(p.Day,nGrp*nB,1), ...
    nW(:), nA(:), mW(:), mA(:), stat(:), df(:), praw(:), pcorr(:), ...
    repmat(string(p.Correction),nGrp*nB,1), repmat(m,nGrp*nB,1), test(:), ...
    'VariableNames',{'Variable','Binning','AgeBin','Day','n_WT','n_APP', ...
    'mean_WT','mean_APP','stat','df','p_raw','p_corr','correction','m','test'});

if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    saveas(f,fullfile(p.SaveDir,sprintf('Panel_Day%d_StrategyGroups_%s.png',p.Day,char(binCol))));
    close(f);
end
end