function PLOT_Tg_GroupUsage_Panel(T, binCol, varargin)
% PLOT_Tg_GroupUsage_Panel(T, 'AgeGroup', 'SaveDir', dir_out)
%
% One slide-ready figure: tight panels, one per genotype (WT, APP). Each panel
% is a grouped bar of the three strategy groups across the age bins of binCol,
% for Day-`Day` usage. Descriptive (the WT-vs-APP test lives in the Day-4
% strategy-group panel).
%
% varargin: 'Day'(4), 'Genotypes'(["WT","APP"]), 'Groups', 'GroupLabels',
%           'Colors'(group), 'YLabel', 'SaveDir'.
% SS 2026

p.Day=4; p.Genotypes=["WT","APP"];
p.Groups={'PlatformIndependent','Procedural','Allocentric'};
p.GroupLabels={'Platform-Independent','Procedural','Allocentric'};
p.Colors=[0.3961 0.2627 0.1294; 1.0000 0.7020 0.4000; 0 0.5020 0];
p.YLabel='Probability for strategy group'; p.SaveDir='';
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end

req=[{'Animal','Genotype','Day',char(binCol)}, p.Groups];
have=T.Properties.VariableNames; miss=req(~ismember(req,have));
if ~isempty(miss)
    error('PLOT_Tg_GroupUsage_Panel:missingCols','T missing column(s): %s',strjoin(miss,', '));
end

G=string(p.Genotypes); nG=numel(G); nGrp=numel(p.Groups);
Td=T(T.Day==p.Day,:);

f=figure('Position',[60 80 1100 480]);
tl=tiledlayout(1,nG,'TileSpacing','tight','Padding','tight');
hKeep=gobjects(nGrp,1);

for gi=1:nG                                      % LOOP: one tile per genotype
    Tgi=Td(string(Td.Genotype)==G(gi),:);
    rat=groupsummary(Tgi,{char(binCol),'Animal'},'mean',p.Groups);
    bins=categories(removecats(Tgi.(binCol)));
    nB=numel(bins);
    vals=cell(nB,nGrp);
    for b=1:nB
        inBin=rat.(binCol)==bins{b};
        for s=1:nGrp
            v=rat.(['mean_' p.Groups{s}])(inBin);
            vals{b,s}=v(~isnan(v));
        end
    end
    ax=nexttile;
    hBar=HELP_DrawGroupedBars(ax,vals,p.Colors,string(p.GroupLabels),bins,0.45,0.10);
    if gi==1, hKeep=hBar; end
    title(ax,char(G(gi)),'FontWeight','bold');
    ylim(ax,[0 1.05]);
    if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(12,12); end
end

xlabel(tl,'Age bin','FontSize',14,'FontWeight','bold');
ylabel(tl,p.YLabel,'FontSize',14,'FontWeight','bold');
title(tl,sprintf('Strategy-group usage  -  Day %d',p.Day),'FontSize',16,'FontWeight','bold');
legend(hKeep,p.GroupLabels,'Orientation','horizontal','Location','southoutside');

if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    saveas(f,fullfile(p.SaveDir,sprintf('Panel_Day%d_GroupUsage_%s.png',p.Day,char(binCol))));
    close(f);
end
end