classdef SymmetricDirichlet
    properties
        conc           % Concentration parameter
        dim            % Dimension
        log_normalizer % Log normalizer
    end

    methods
        % Constructor
        function obj = SymmetricDirichlet(conc, dim)
            obj.conc = conc;
            obj.dim = dim;
            obj.log_normalizer = dim * gammaln(conc) - gammaln(dim * conc);
        end

        % Random number generation
        function x = rand(obj, dim)
            % Updates the input vector x with random values from a Symmetric Dirichlet distribution

            % Check if the size of x matches the dimension of the distribution
            if nargin == 1
                dim = obj.dim;
            end

            if isscalar(dim)
                dim = [1,dim];
            end

            % Generate samples from a Gamma distribution
            x = gamrnd(obj.conc, 1, dim);

            % Normalize the samples
            x = x ./ sum(sum(x));
        end

        % Posterior calculation
        function post = posterior(obj, counts)
            % Assuming counts is a vector of the same length as obj.dim
            assert(length(counts) == obj.dim, 'Counts vector must match the dimension of the Dirichlet distribution.');
            post = Dirichlet(counts + obj.conc); % Dirichlet class needs to be defined
        end

        % Log probability density function
        function val = logpdf(obj, p)
            assert(length(p) == obj.dim, 'Probability vector must match the dimension of the Dirichlet distribution.');
            val = sum(log(p) * (obj.conc - 1)) - obj.log_normalizer;
        end

        % Length method
        function l = length(obj)
            l = obj.dim;
        end
    end
end
