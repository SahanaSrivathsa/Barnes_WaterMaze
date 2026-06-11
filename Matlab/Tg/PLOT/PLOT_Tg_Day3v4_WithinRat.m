function S = PLOT_Tg_Day3v4_WithinRat(T, varName, varargin)
% S = PLOT_Tg_Day3v4_WithinRat(T, 'Procedural', 'GroupBy','AgeGroup', ...
%                              'Levels',["Young","Mid","Old"], 'Colors',greyShades)
%
% Single-figure within-rat Day d1->d2 plot, coloured/averaged by GroupBy.
% Returns a paired-test stats table (one row per level).
%
% INPUT
%   T, varName : Build_Tg_StrategyTable output + variable to plot.
%   varargin (name-value):
%     'GroupBy'   factor column (default 'Genotype')
%     'Levels'    string array level order (default observed)
%     'Colors'    nLevel x 3 (default grey ramp light->dark)
%     'Filter'    {col,value} subset before plotting (default {})
%     'Days'      [d1 d2] (default [3 4])
%     'IndivAlpha' default 0.25
%     'Test'      'ttest' (paired, default) | 'signrank'
%     'YLabel','TitleTag','SaveDir'
% OUTPUT
%   S : table {Variable,GroupBy,Tag,Level,n,dayA,dayB,mean_dayA,mean_dayB,stat,df,p,test}
%
% SS 2026


p.GroupBy='Genotype'; p.Levels=strings(0); p.Colors=[]; p.Filter={};
p.Days=[3 4]; p.IndivAlpha=0.25; p.Test='ttest';
p.YLabel=char(varName); p.TitleTag=''; p.SaveDir='';
for k=1:2:numel(varargin), p.(varargin{k})=varargin{k+1}; end
 
lab = HELP_StrategyLabel(varName); tok = strrep(lab,' ','_');
 
req={'Animal','Day',p.GroupBy,char(varName)};
if ~isempty(p.Filter), req=[req,p.Filter(1)]; end
have=T.Properties.VariableNames;
miss=req(~ismember(req,have));
if ~isempty(miss)
    error('PLOT_Tg_Day3v4_WithinRat:missingCols','T missing column(s): %s',strjoin(miss,', '));
end
 
if ~isempty(p.Filter), T=T(string(T.(p.Filter{1}))==string(p.Filter{2}),:); end
d1=p.Days(1); d2=p.Days(2);
T=T(ismember(T.Day,[d1 d2]),:);
S=table();
if isempty(T)
    warning('PLOT_Tg_Day3v4_WithinRat:noData','No rows after filtering for %s.',char(varName));
    return
end
 
grpStr=string(T.(p.GroupBy));
if isempty(p.Levels), L=unique(grpStr); L=L(~ismissing(L)&L~=""); else, L=string(p.Levels); end
nL=numel(L);
if isempty(p.Colors)
    if nL==1, p.Colors=[0.4 0.4 0.4]; else, g=linspace(0.72,0.12,nL)'; p.Colors=[g g g]; end
end
 
rat=groupsummary(T,{'Animal',p.GroupBy,'Day'},'mean',char(varName));
ratCol=['mean_' char(varName)];
 
f=figure('Position',[100 100 560 620]);
[hMean,pairByLevel]=HELP_DrawDay3v4Axes(gca,rat,ratCol,p.GroupBy,L,[d1 d2],p.Colors,p.IndivAlpha);
 
if isempty(p.TitleTag), tag=char(p.GroupBy); else, tag=char(p.TitleTag); end
xlabel('Day'); ylabel(p.YLabel);
title(sprintf('%s  -  %s  (Day %d vs %d)',lab,tag,d1,d2),'Interpreter','none');
legend(hMean,cellstr(L),'Location','northeastoutside');
if exist('pubify_figure_axis_robust','file'), pubify_figure_axis_robust(14,14); end
 
% ---- paired stats per level 
Variable=repmat(string(lab),nL,1);
GroupBy =repmat(string(p.GroupBy),nL,1);
Tag     =repmat(string(tag),nL,1);
Level   =string(L(:));
dayA=repmat(d1,nL,1); dayB=repmat(d2,nL,1);
n=zeros(nL,1); mA=nan(nL,1); mB=nan(nL,1); stat=nan(nL,1); df=nan(nL,1); pv=nan(nL,1); test=strings(nL,1);
for li=1:nL                                      % LOOP: per level (tiny)
    P=pairByLevel{li}; P=P(all(~isnan(P),2),:); n(li)=size(P,1);
    if n(li)>=1, mA(li)=mean(P(:,1)); mB(li)=mean(P(:,2)); end
    if n(li)>=2
        if strcmpi(p.Test,'signrank')
            [pp,~,st]=signrank(P(:,1),P(:,2)); pv(li)=pp; test(li)="signrank";
            if isfield(st,'zval'), stat(li)=st.zval; end
        else
            [~,pp,~,st]=ttest(P(:,1),P(:,2)); pv(li)=pp; stat(li)=st.tstat; df(li)=st.df; test(li)="paired_ttest";
        end
    end
end
S=table(Variable,GroupBy,Tag,Level,n,dayA,dayB,mA,mB,stat,df,pv,test, ...
    'VariableNames',{'Variable','GroupBy','Tag','Level','n','dayA','dayB','mean_dayA','mean_dayB','stat','df','p','test'});
 
if ~isempty(p.SaveDir)
    if ~exist(p.SaveDir,'dir'), mkdir(p.SaveDir); end
    saveas(f,fullfile(p.SaveDir,sprintf('Day%dv%d_%s_%s.png',d1,d2,tok,tag)));
    close(f);
end
end