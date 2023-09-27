function [Y, AUC] = getAUC(V,idx1,idx2,AUCPrecision)
    % V: vector with values
    % idx1: indices of class 1
    % idx2: indices of class 2
    % AUCPrecision: precision of the ROC

    labels = [zeros(1, numel(idx1)), ones(1, numel(idx2))];
    scores = [V(idx1)', V(idx2)'];
    [X, Y, ~, AUC] = perfcurve(labels, scores, 1,'XVals',AUCPrecision);

    % find closest precision value for every X
    [~,idx] = min(abs(AUCPrecision-X),[],1);
    Y = Y(idx);

end