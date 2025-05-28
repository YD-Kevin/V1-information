function J = computeLinearFI(tuningCurve, responses, stimulus)

dFds = zeros(size(tuningCurve));
for i = 1:size(tuningCurve, 1)
    dFds(i,:) = gradient(tuningCurve(i,:), stimulus);
end


covMatrix = cov(responses');  


epsilon = 1e-5;
covMatrix = covMatrix + epsilon * eye(size(covMatrix));


J = zeros(1, length(stimulus));
for k = 1:length(stimulus)
    df = dFds(:, k);
    J(k) = df' * (covMatrix \ df);  % J = (df/ds)^T Σ^{-1} (df/ds)
end



end
