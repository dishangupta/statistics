% This functions samples from the log-linear model
% log p(x) = beta0 + 3(x1x2 + x2x3 + x3x4 + x4x5)
% and fits a maximum likelihood estimator on 
% log p(x) = beta0 + \sum_j(betaj*xj) + \sum_{k<l}(betakl*xk*xl)

function [b, dev, stats] = SMLHW4Q5(N)
    % ----- Generate probability distribution -----

    % Initialize probability and domain
    D = [0:2^5-1]';
    X = rem(floor(D*pow2(-(5-1):0)),2);
    P = zeros(size(X, 1), 1);

    syms x1 x2 x3 x4 x5;
    px = exp(5 + (x1*x2 + x2*x3 + x3*x4 + x4*x5)); % Choose beta0 to be some arbitrary value
    for j = 1:size(X, 1)
        p = subs(px, {x1, x2, x3, x4, x5}, X(j, :));
        P(j) = double(p);
    end

    % Normalize P 
    Pdash = P;
    P = P/sum(P);

    % Sample N random vectors from P
    S = mnrnd(N, P);
    
    % Fit MLE for model: 
    % log p(x) = beta0 + \sum_j(betaj*xj) + \sum_{k<l}(betakl*xk*xl)
    
    Xtrain = [];
    yTrain = [];
    
    numSamples = 0;
    for j = 1:size(S, 2)
        numSamples = S(j);

        while (numSamples ~= 0)
            f12 = X(j, 1)*X(j, 2);
            f13 = X(j, 1)*X(j, 3);
            f14 = X(j, 1)*X(j, 4);
            f15 = X(j, 1)*X(j, 5);
            f23 = X(j, 2)*X(j, 3);
            f24 = X(j, 2)*X(j, 4);
            f25 = X(j, 2)*X(j, 5);
            f34 = X(j, 3)*X(j, 4);
            f35 = X(j, 3)*X(j, 5);
            f45 = X(j, 4)*X(j, 5);
            Xtrain = [Xtrain; X(j, :), f12, f13, f14, f15, f23, f24, f25, f34, f35, f45];
            yTrain = [yTrain; Pdash(j)]; % To change to P
            numSamples = numSamples - 1;
        end

    end
    
    [b,dev,stats] = glmfit(Xtrain, yTrain, 'poisson', 'link', 'log');
     
    
    
end

