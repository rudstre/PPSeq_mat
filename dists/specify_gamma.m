function rg = specify_gamma(meanVal, varVal)
    % Calculates Gamma distribution parameters from mean and variance
    % and returns a RateGamma structure

    beta = meanVal / varVal;
    alpha = meanVal * beta;

    rg = RateGamma(alpha, beta);
end