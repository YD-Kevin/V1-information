function bestT=findSuffStats(data, candidateFuncs, maxk)

n=length(data);
numCandidates = length(candidateFuncs);
aicValues = [];
selectedIndices = {};

for k=1:maxk
    combinations=nchoosek(1:numCandidates, k);
    for i = 1:size(combinations, 1)
        indices=combinations(i, :);
        T_funcs=candidateFuncs(indices);
        

        TData=zeros(n, length(indices));
        for j=1:n
            TRow=[];
            for f=T_funcs
                TRow=[TRow, f{1}(data(j))];
            end
            TData(j, :) = TRow;
        end
        

        muT = mean(TData)';
        sigmaT = cov(TData);
        

        logLikelihood = -0.5 * n * log(det(sigmaT)) - 0.5 * trace((TData - muT') / sigmaT * (TData - muT')');
        num_params = k + k*(k+1)/2; 
        aic = 2*num_params - 2*logLikelihood;
        
        aicValues = [aicValues; aic];
        selectedIndices{end+1} = indices;
    end
end

[~, idx] = min(aicValues);
bestT = selectedIndices{idx};
end

