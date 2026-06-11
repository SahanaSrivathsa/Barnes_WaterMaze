function T = Build_Tg_StrategyTable(rawT, varargin)
% Reshape  RTrack strategy output to long form
% Adds cols Animal, Day, Sex, Genotype, AgeGroup  Young/Mid/Old)
% AgeMonth
% Groups (Sex x Genotype),
% the per-trial strategy-group sums (entropy),
% SS 2026
%% PARAMS
p.StrategyVars = {'thigmotaxis','circling','randomPath','scanning', ...
'chaining','directedSearch','correctedPath','directPath'};
p.GroupDef     = { ...
'PlatformIndependent', {'thigmotaxis','circling','randomPath'}; ...
'Procedural',          {'scanning','chaining'}; ...
'Allocentric',         {'directedSearch','correctedPath','directPath'} };
p.AgeEdges    = [3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 17.5 22.5];
p.AgeLabels   = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12mo","13mo","14mo","15-16mo","20-21mo"];
p.YoungMidCut = 9;
p.MidOldCut   = 15;
for k = 1:2:numel(varargin)
    p.(varargin{k}) = varargin{k+1};
end
% Validate cols
req  = [{'x_TargetID','x_Day','Age','Sex','Genotype'}, p.StrategyVars];
have = rawT.Properties.VariableNames;
miss = req(~ismember(req, have));
if ~isempty(miss)
    error('Build_Tg_StrategyTable:missingCols', ...
'rawT missing required column(s): %s', strjoin(miss, ', '));
end
n = height(rawT);
T = table;
T.Animal = string(rawT.x_TargetID);
T.Day    = double(rawT.x_Day);
% Age in months
ageRaw = rawT.Age;
if iscell(ageRaw)
    ageNum = cellfun(@(x) str2double(regexp(num2str(x),'[\d.]+','match','once')), ageRaw);
elseif isstring(ageRaw) || ischar(ageRaw)
    ageNum = arrayfun(@(x) str2double(regexp(x,'[\d.]+','match','once')), string(ageRaw));
else
    ageNum = double(ageRaw);
end
% Sex + Genotype  (genotype lives in the Genotype column: APP/+ , WT)
T.Sex = string(rawT.Sex);
genoIn   = string(rawT.Genotype);
genotype = strings(n,1);
genotype(contains(genoIn,"APP",'IgnoreCase',true)) = "APP";
genotype(strcmpi(genoIn,"WT"))                      = "WT";
T.Genotype = genotype;
% Broad age
coarse = repmat("Old", n, 1);
coarse(ageNum <= p.YoungMidCut)                                 = "Young";
coarse(ageNum >  p.YoungMidCut & ageNum <= p.MidOldCut)         = "Mid";
coarse(isnan(ageNum))                                           = missing;
T.AgeGroup = categorical(coarse, ["Young","Mid","Old"], 'Ordinal', true);
%Fine month bins
T.AgeMonth = AssignAgeBin(ageNum, p.AgeEdges, p.AgeLabels);
% Groupss (Sex x Genotype)
T.Group4 = categorical(T.Sex + " " + T.Genotype, ["F WT","F APP","M WT","M APP"]);
% Prob cols
for s = 1:numel(p.StrategyVars)
    T.(p.StrategyVars{s}) = double(rawT.(p.StrategyVars{s}));
end
% ---- strategy-group sums ----
nGrp     = size(p.GroupDef,1);
grpNames = string(p.GroupDef(:,1));
grpSums  = zeros(n, nGrp);
for g = 1:nGrp                            % LOOP: per group; member sets differ in length
    members = p.GroupDef{g,2};
    bad     = members(~ismember(members, p.StrategyVars));
if ~isempty(bad)
        error('Build_Tg_StrategyTable:badMember', ...
'GroupDef "%s" references unknown strategy/strategies: %s', ...
            grpNames(g), strjoin(bad, ', '));
end
    grpSums(:,g)        = sum(rawT{:, members}, 2);
    T.(char(grpNames(g))) = grpSums(:,g);
end
% ---- per-trial Shannon entropies (bits) ----
pAll  = T{:, p.StrategyVars};
T.Entropy_All   = -sum(pAll  .* log2(pAll  + eps), 2);
T.Entropy_Group = -sum(grpSums .* log2(grpSums + eps), 2);
end