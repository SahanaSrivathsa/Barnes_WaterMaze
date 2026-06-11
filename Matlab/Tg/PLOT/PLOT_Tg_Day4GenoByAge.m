function S = PLOT_Tg_Day4GenoByAge(T, varName, binCol, varargin)
% S = PLOT_Tg_Day4GenoByAge(T, 'thigmotaxis', 'AgeGroup', 'SaveDir', dir_out)
%
% Single-figure grouped bar: WT vs APP for one variable on one day, age bins of
% binCol on the x-axis (rat = unit, mean +/- SEM, rat dots). WT-vs-APP test per
% bin marked with sigstar. Returns a stats table (one row per age bin).
%
% varargin: 'Day'(4), 'Genotypes'(["WT","APP"]), 'Colors', 'Test'('ttest2'|'ranksum'),
%           'YLabel', 'SaveDir'.
% SS 2026

p.Day=4; p.Genotypes=["WT","APP"];
p.Colors=[0.50 0.70 1.00; 0.00 0.20 0.70];
p.Test='ttest2'; p.YLabel=char(varName); p.SaveDir='';
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end

lab = HELP_StrategyLabel(varName); tok = strrep(lab,' ','_');

req={'Animal','Day','Genotype',char(binCol),char(varName)};
have=T.Properties.VariableNames; miss=req(~ismember(req,have));
if ~isempty(miss)
    error('PLOT_Tg_Day4GenoByAge:missingCols','T missing column(s): %s',strjoin(miss,', '));
end

Td  = T(T.Day==p.Day,:);
rat = groupsummary(Td,{char(binCol),'Genotype','Animal'},'mean',char(varName));
ratCol=['mean_' char(varName)];
bins=categories(removecats(Td.(binCol)));
G=string(p.Genotypes); nB=numel(bins); nG=numel(G);

vals=cell(nB,nG);
for b=1:nB                                       % LOOP: bins x genotypes (small)
    for g=1:nG
        v=rat.(ratCol)(rat.(binCol)==bins{b} & string(rat.Genotype)==G(g));
        vals{b,g}=v(~isnan(v));
    end
end

f=figure('Position',[100 100 1000 560]);
[hBar,ctr]=HELP_DrawGroupedBars(gca,vals,p.Colors,G,bins,0.75,0.12);

% WT vs APP per bin -> stats + sigstar
Variable=repmat(string(lab),nB,1); Binning=repmat(string(binCol),nB,1);
AgeBin=string(bins(:)); Day=repmat(p.Day,nB,1);
nW=zeros(nB,1); nA=zeros(nB,1); mW=nan(nB,1); mA=nan(nB,1);
stat=nan(nB,1); df=nan(nB,1); pv=nan(nB,1); test=strings(nB,1);
sigPairs={}; sigP=[];
for b=1:nB
    a=vals{b,1}; c=vals{b,2}; nW(b)=numel(a); nA(b)=numel(c);
    if ~isempty(a), mW(b)=mean(a); end
    if ~isempty(c), mA(b)=mean(c); end
    if numel(a)>1 && numel(c)>1
        if strcmpi(p.Test,'ranksum')
            [pp,~,st]=ranksum(a,c); test(b)="ranksum"; if isfield(st,'zval'), stat(b)=st.zval; end
        else
            [~,pp,~,st]=ttest2(a,c,'Vartype','unequal'); stat(b)=st.tstat; df(b)=st.df; test(b)="ttest2_Welch";
        end
        pv(b)=pp;
        if pp<0.05, sigPairs{end+1}=[ctr(b,1),ctr(b,2)]; sigP(end+1)=pp; end %#ok<AGROW>
    end
end
if ~isempty(sigPairs) && exist('sigstar','file'), sigstar(sigPairs,sigP); end

xlabel('Age bin'); ylabel(p.YLabel);
title(sprintf('%s  -  Day %d  (WT vs APP)',lab,p.Day),'Interpreter','none');
legend(hBar,cellstr(G),'Location','northeastoutside');
if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(14,14); end
hold off;

S=table(Variable,Binning,AgeBin,Day,nW,nA,mW,mA,stat,df,pv,test, ...
    'VariableNames',{'Variable','Binning','AgeBin','Day','n_WT','n_APP','mean_WT','mean_APP','stat','df','p','test'});

if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    saveas(f,fullfile(p.SaveDir,sprintf('MWM_Day%d_%s_%s.png',p.Day,tok,char(binCol))));
    close(f);
end
end