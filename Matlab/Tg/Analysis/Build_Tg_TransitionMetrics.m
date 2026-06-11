function T = Build_Tg_TransitionMetrics(rawT, varargin)
% Build per-rat-per-day within-day strategy transition matrices matrices
% make foll metrics
%   TransEntropy : Shannon entropy via matrixEntropy 
%   SpectralGap  :  the gap between the largest and secondlargest
%  eigenvalue magnitudes of the transition matrix.
%
% SS 2026

% Input params
p.StrategyVars   = {'thigmotaxis','circling','randomPath','scanning', ...
                    'chaining','directedSearch','correctedPath','directPath'};
p.StateSource    = 'argmax';
p.StateColumn    = 'strategy';
p.AgeEdges       = [3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 17.5 22.5];
p.AgeLabels      = ["4mo","5mo","6mo","7mo","8mo","9mo","10mo","11mo","12mo","13mo","14mo","15-16mo","20-21mo"];
p.YoungMidCut    = 9;
p.MidOldCut      = 15;
p.MinTrials      = 2;
p.MakeStochastic = false;
for k = 1:2:numel(varargin)              
    p.(varargin{k}) = varargin{k+1};
end

% Validate cols
req = {'x_TargetID','x_Day','x_Trial','Age','Sex','Genotype'};
if strcmpi(p.StateSource,'argmax')
    req = [req, p.StrategyVars];
else
    req = [req, {p.StateColumn}];
end
have = rawT.Properties.VariableNames;
miss = req(~ismember(req, have));
if ~isempty(miss)
    error('Build_Tg_TransitionMetrics:missingCols', ...
        'rawT missing required column(s): %s', strjoin(miss, ', '));
end

K = numel(p.StrategyVars);

% Age
ageRaw = rawT.Age;
if iscell(ageRaw)
    ageNum = cellfun(@(x) str2double(regexp(num2str(x),'[\d.]+','match','once')), ageRaw);
elseif isstring(ageRaw) || ischar(ageRaw)
    ageNum = arrayfun(@(x) str2double(regexp(x,'[\d.]+','match','once')), string(ageRaw));
else
    ageNum = double(ageRaw);
end

% Discrete state per trial
if strcmpi(p.StateSource,'argmax')
    [~, stateVec] = max(rawT{:, p.StrategyVars}, [], 2);   
else
    sc = rawT.(p.StateColumn);
    if ~isnumeric(sc), sc = double(categorical(sc)); end      
    stateVec = double(sc);
    if any(stateVec < 1 | stateVec > K | mod(stateVec,1)~=0)
        error('Build_Tg_TransitionMetrics:badState', ...
            'StateColumn "%s" must be integers in 1..%d.', p.StateColumn, K);
    end
end

animals = string(rawT.x_TargetID);
days    = double(rawT.x_Day);
trials  = double(rawT.x_Trial);
sexAll  = string(rawT.Sex);
genoAll = string(rawT.Genotype);

uAnimals = unique(animals, 'stable');
uDays    = unique(days);
nA = numel(uAnimals);
nD = numel(uDays);

% preallocate
maxRows   = nA * nD;
outAnimal = strings(maxRows,1);
outDay    = nan(maxRows,1);
outSex    = strings(maxRows,1);
outGeno   = strings(maxRows,1);
outAge    = nan(maxRows,1);
outEnt    = nan(maxRows,1);
outGap    = nan(maxRows,1);
cnt = 0;

for ai = 1:nA                         
    am   = animals == uAnimals(ai);
    aRow = find(am, 1);
    % per-animal 
    aSex     = sexAll(aRow);
    aGenoRaw = genoAll(aRow);
    if contains(aGenoRaw, "APP", 'IgnoreCase', true)
        aGeno = "APP";
    elseif strcmpi(aGenoRaw, "WT")
        aGeno = "WT";
    else
        aGeno = missing;
    end
    aAge = ageNum(aRow);

    for di = 1:nD                        
        dm = am & (days == uDays(di));
        if ~any(dm), continue; end
        cnt = cnt + 1;
        outAnimal(cnt) = uAnimals(ai);
        outDay(cnt)    = uDays(di);
        outSex(cnt)    = aSex;
        outGeno(cnt)   = aGeno;
        outAge(cnt)    = aAge;

        % ordered within-day state sequence
        [~, ord] = sort(trials(dm));
        st       = stateVec(dm);
        st       = st(ord);
        if numel(st) < p.MinTrials, continue; end

        % consecutive-trial transition
        cm = zeros(K);
        for t = 1:numel(st)-1          
            cm(st(t), st(t+1)) = cm(st(t), st(t+1)) + 1;
        end
        if ~any(cm(:)), continue; end

        % row normalize
        rs = sum(cm, 2);
        nz = rs > 0;
        cm(nz,:) = cm(nz,:) ./ rs(nz);
        if p.MakeStochastic && any(~nz)
            cm(~nz,:) = 1 / K;           
        end

        
        outEnt(cnt) = matrixEntropy(cm);
        ev  = eig(cm);
        lam = sort(abs(ev), 'descend');
        if numel(lam) >= 2, outGap(cnt) = 1 - lam(2); end
    end
end

% trim
idx = 1:cnt;
T = table;
T.Animal       = outAnimal(idx);
T.Day          = outDay(idx);
T.Sex          = outSex(idx);
T.Genotype     = outGeno(idx);
T.TransEntropy = outEnt(idx);
T.SpectralGap  = outGap(idx);


coarse = repmat("Old", cnt, 1);
coarse(outAge(idx) <= p.YoungMidCut)                                 = "Young";
coarse(outAge(idx) >  p.YoungMidCut & outAge(idx) <= p.MidOldCut)    = "Mid";
coarse(isnan(outAge(idx)))                                           = missing;
T.AgeGroup = categorical(coarse, ["Young","Mid","Old"], 'Ordinal', true);
T.AgeMonth = AssignAgeBin(outAge(idx), p.AgeEdges, p.AgeLabels);
T.Group4   = categorical(T.Sex + " " + T.Genotype, ["F WT","F APP","M WT","M APP"]);
end