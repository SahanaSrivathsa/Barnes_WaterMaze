function S = PLOT_Tg_Day3v4_AgeSubplots(T, varName, varargin)

% SS 2026


p.AgeCol='AgeGroup'; p.AgeLevels=["Young","Mid","Old"]; p.Genotypes=["WT","APP"];
p.Colors=[0.50 0.70 1.00; 0.00 0.20 0.70]; p.Days=[3 4]; p.IndivAlpha=0.25;
p.Test='ttest'; p.YLabel=char(varName); p.SaveDir='';
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end
 
lab = HELP_StrategyLabel(varName); tok = strrep(lab,' ','_');
 
req={'Animal','Day','Genotype',p.AgeCol,char(varName)};
have=T.Properties.VariableNames;
miss=req(~ismember(req,have));
if ~isempty(miss)
    error('PLOT_Tg_Day3v4_AgeSubplots:missingCols','T missing column(s): %s',strjoin(miss,', '));
end
 
d1=p.Days(1); d2=p.Days(2);
G=string(p.Genotypes); nG=numel(G);
ages=string(p.AgeLevels);
present=ages(ismember(ages, unique(string(T.(p.AgeCol)))));   % keep populated, in order
nA=numel(present);
if nA==0
    warning('PLOT_Tg_Day3v4_AgeSubplots:noData','No populated age levels for %s.',char(varName));
    S=table(); return
end
 
% preallocate stats
nRows=nA*nG;
Variable=strings(nRows,1); AgeGroup=strings(nRows,1); Genotype=strings(nRows,1);
n=zeros(nRows,1); dayA=repmat(d1,nRows,1); dayB=repmat(d2,nRows,1);
mA=nan(nRows,1); mB=nan(nRows,1); stat=nan(nRows,1); df=nan(nRows,1); pv=nan(nRows,1); test=strings(nRows,1);
kk=0;
 
f=figure('Position',[60 80 360*nA 460]);
tl=tiledlayout(1,nA,'TileSpacing','tight','Padding','tight');
hMeanKeep=gobjects(nG,1);
 
for a=1:nA                                       % LOOP: one tile per age group
    ax=nexttile;
    Ta=T(string(T.(p.AgeCol))==present(a) & ismember(T.Day,[d1 d2]),:);
    rat=groupsummary(Ta,{'Animal','Genotype','Day'},'mean',char(varName));
    ratCol=['mean_' char(varName)];
    [hMean,pairByLevel]=HELP_DrawDay3v4Axes(ax,rat,ratCol,'Genotype',G,[d1 d2],p.Colors,p.IndivAlpha);
    title(ax,char(present(a)),'FontWeight','bold');
    if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(12,12); end
    if a==1, hMeanKeep=hMean; end
 
    for gi=1:nG                                  % LOOP: per genotype (stats)
        P=pairByLevel{gi}; P=P(all(~isnan(P),2),:);
        kk=kk+1;
        Variable(kk)=string(lab); AgeGroup(kk)=present(a); Genotype(kk)=G(gi);
        n(kk)=size(P,1);
        if n(kk)>=1, mA(kk)=mean(P(:,1)); mB(kk)=mean(P(:,2)); end
        if n(kk)>=2
            if strcmpi(p.Test,'signrank')
                [pp,~,st]=signrank(P(:,1),P(:,2)); pv(kk)=pp; test(kk)="signrank";
                if isfield(st,'zval'), stat(kk)=st.zval; end
            else
                [~,pp,~,st]=ttest(P(:,1),P(:,2)); pv(kk)=pp; stat(kk)=st.tstat; df(kk)=st.df; test(kk)="paired_ttest";
            end
        end
    end
end
 
xlabel(tl,'Day','FontSize',14,'FontWeight','bold');
ylabel(tl,p.YLabel,'FontSize',14,'FontWeight','bold');
title(tl,sprintf('%s  (Day %d vs %d)',lab,d1,d2),'FontSize',16,'FontWeight','bold','Interpreter','none');
legend(hMeanKeep,cellstr(G),'Orientation','horizontal','Location','southoutside');
 
S=table(Variable,AgeGroup,Genotype,n,dayA,dayB,mA,mB,stat,df,pv,test, ...
    'VariableNames',{'Variable','AgeGroup','Genotype','n','dayA','dayB','mean_dayA','mean_dayB','stat','df','p','test'});
 
if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    saveas(f,fullfile(p.SaveDir,sprintf('Day%dv%d_%s_AgeSubplots.png',d1,d2,tok)));
    close(f);
end
end
 