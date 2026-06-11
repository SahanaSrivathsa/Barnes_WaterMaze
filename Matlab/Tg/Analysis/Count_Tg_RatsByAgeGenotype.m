function C = Count_Tg_RatsByAgeGenotype(T)

% Count UNIQUE rats per age group per genotype, for both the coarse
% (AgeGroup: Young/Mid/Old) and fine (AgeMonth: month bins) binnings.
% One rat = one unique Animal id (each maps to a single age/genotype).
%
%

req  = {'Animal','Genotype','AgeGroup','AgeMonth'};
have = T.Properties.VariableNames;
miss = req(~ismember(req, have));
if ~isempty(miss)
    error('Count_Tg_RatsByAgeGenotype:missingCols', ...
        'T missing column(s): %s', strjoin(miss, ', '));
end

% one row per rat
[~, ia] = unique(string(T.Animal), 'stable');
R    = T(ia, :);
geno = string(R.Genotype);

gset = unique(geno);
gset = gset(~ismissing(gset) & gset ~= "");
nG   = numel(gset);

agesC = categories(removecats(R.AgeGroup));   % present coarse bins, ordinal order
agesF = categories(removecats(R.AgeMonth));   % present month bins, ordinal order

% preallocate
nRows    = (numel(agesC) + numel(agesF)) * nG;
Binning  = strings(nRows,1);
AgeBin   = strings(nRows,1);
Genotype = strings(nRows,1);
nRats    = zeros(nRows,1);

binDefs = {'Coarse', string(R.AgeGroup), agesC; ...
           'Fine',   string(R.AgeMonth), agesF};
k = 0;
for b = 1:size(binDefs,1)                    
    binName = binDefs{b,1};
    colVals = binDefs{b,2};
    ages    = binDefs{b,3};
    for ai = 1:numel(ages)
        for gi = 1:nG
            k = k + 1;
            Binning(k)  = binName;
            AgeBin(k)   = string(ages{ai});
            Genotype(k) = gset(gi);
            nRats(k)    = sum(colVals == string(ages{ai}) & geno == gset(gi));
        end
    end
end

C = table(Binning, AgeBin, Genotype, nRats);
end