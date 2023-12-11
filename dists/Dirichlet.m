classdef Dirichlet
    properties
        alpha   % Parameter vector of the Dirichlet distribution
    end

    methods
        % Constructor
        function obj = Dirichlet(alpha)
            obj.alpha = alpha;
        end

        % Method to generate random samples from Dirichlet distribution
        function samples = rand(obj)
            samples = gamrnd(obj.alpha, 1);
            samples = samples ./ sum(samples);
        end

        % Log probability density function
        function val = logpdf(obj, x)
            % Validate the input
            assert(all(x >= 0) && abs(sum(x) - 1) < 1e-10, 'Input must be a probability vector.');
            
            % Compute the log PDF
            val = sum((obj.alpha - 1) .* log(x)) + gammaln(sum(obj.alpha)) - sum(gammaln(obj.alpha));
        end

        % Method to compute the mean of the distribution
        function mean_val = mean(obj)
            mean_val = obj.alpha / sum(obj.alpha);
        end

        % Variance of the distribution (returns a matrix)
        function var_matrix = variance(obj)
            alpha0 = sum(obj.alpha);
            var_matrix = zeros(length(obj.alpha));
            for i = 1:length(obj.alpha)
                for j = 1:length(obj.alpha)
                    if i == j
                        var_matrix(i, j) = obj.alpha(i) * (alpha0 - obj.alpha(i)) / (alpha0^2 * (alpha0 + 1));
                    else
                        var_matrix(i, j) = -obj.alpha(i) * obj.alpha(j) / (alpha0^2 * (alpha0 + 1));
                    end
                end
            end
        end
    end
end
