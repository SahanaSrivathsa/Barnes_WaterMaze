function PLOT_Tg_GroupUsageByAge(T, binCol, varargin)
% Grouped bar of strategy-GROUP usage: age bins of binCol on the x-axis, the
% three strategy groups as the bar series, for ONE genotype. If 'Day' is set,
% only that day is used (per-rat day mean); otherwise usage is collapsed over
% all days. Across-rat mean +/- SEM with rat dots overlaid.
%
% INPUT
%   T      : Build_Tg_StrategyTable output (needs the group-sum columns).
%   binCol : 'AgeGroup' or 'AgeMonth'.
%   
% SS 2026

p.Genotype="WT"; p.Day=[];
p.Groups={'PlatformIndependent','Procedural','Allocentric'};
p.GroupLabels={'Platform-Independent','Procedural','Allocentric'};
p.Colors=[0.3961 0.2627 0.1294; 1.0000 0.7020 0.4000; 0 0.5020 0];
p.YLabel='Probability for strategy group'; p.SaveDir='';
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end
 
req=[{'Animal','Genotype','Day',char(binCol)}, p.Groups];
have=T.Properties.VariableNames; miss=req(~ismember(req,have));
if ~isempty(miss)
    error('PLOT_Tg_GroupUsageByAge:missingCols','T missing column(s): %s',strjoin(miss,', '));
end
 
Tg=T(string(T.Genotype)==string(p.Genotype),:);
if ~isempty(p.Day)
    Tg=Tg(Tg.Day==p.Day,:); dayStr=sprintf(' - Day %d',p.Day); fileDay=sprintf('_Day%d',p.Day);
else
    dayStr=''; fileDay='';
end
rat=groupsummary(Tg,{char(binCol),'Animal'},'mean',p.Groups);
bins=categories(removecats(Tg.(binCol)));
nB=numel(bins); nGrp=numel(p.Groups);
 
vals=cell(nB,nGrp);
for b=1:nB                                       % LOOP: bins x groups (small)
    inBin=rat.(binCol)==bins{b};
    for gi=1:nGrp
        v=rat.(['mean_' p.Groups{gi}])(inBin);
        vals{b,gi}=v(~isnan(v));
    end
end
 
f=figure('Position',[100 100 1100 600]);
hBar=HELP_DrawGroupedBars(gca,vals,p.Colors,string(p.GroupLabels),bins,0.45,0.10);
 
xlabel('Age bin'); ylabel(p.YLabel);
title(sprintf('Strategy-group usage  -  %s%s',char(p.Genotype),dayStr),'Interpreter','none');
legend(hBar,p.GroupLabels,'Location','northeastoutside');
ylim([0 1.05]);
if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(14,14); end
hold off;
 
if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    saveas(f,fullfile(p.SaveDir,sprintf('GroupUsage_%s_%s%s.png',char(p.Genotype),char(binCol),fileDay)));
    close(f);
end
end
 