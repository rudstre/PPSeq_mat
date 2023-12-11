classdef RateGamma
    properties
        alpha   % Shape parameter
        beta    % Rate parameter
    end

    methods
        % Constructor
        function obj = RateGamma(alpha, beta)
            obj.alpha = alpha;
            obj.beta = beta;
        end

        % Logarithm of probability density function
        function val = logpdf(obj, x)
            val = (obj.alpha - 1) * log(x) - obj.beta * x + obj.alpha * log(obj.beta) - gammaln(obj.alpha);
        end

        % Generate random sample
        function r = rand(obj)
            r = gamrnd(obj.alpha, 1 / obj.beta);
        end

        % Mean of the distribution
        function m = mean(obj)
            m = obj.alpha / obj.beta;
        end

        % Variance of the distribution
        function v = var(obj)
            v = obj.alpha / (obj.beta^2);
        end
    end

    % Static method for posterior calculation
    methods (Static)
        function post = posterior(count_var, prior)
            post = RateGamma(prior.alpha + count_var, prior.beta + 1);
        end
    end
end
