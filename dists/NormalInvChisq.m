classdef NormalInvChisq
    properties
        k    % Number of pseudo-observations for mean
        m    % Prior mean
        nu   % Degrees of freedom for variance
        s2   % Scale parameter for variance
    end

    methods
        % Constructor
        function obj = NormalInvChisq(k, m, nu, s2)
            obj.k = k;
            obj.m = m;
            obj.nu = nu;
            obj.s2 = s2;
        end

        % Random number generation
        function [mu, sigma2] = rand(obj)
            sigma2 = rand(ScaledInvChisq(obj.nu, obj.s2));
            mu = obj.m + normrnd(0, sqrt(sigma2 / obj.k));
        end

        % Posterior distribution
        function post = posterior(obj, n, sum_x, sumsq_x)
            if n == 0
                post = obj;
                return;
            end
            k_new = obj.k + n;
            nu_new = obj.nu + n;
            m_new = (sum_x / k_new) + (obj.k * obj.m / k_new);
            s2_new = (obj.nu * obj.s2 + (sumsq_x - sum_x * sum_x / n) + (obj.k * n * (sum_x / n - obj.m)^2) / k_new) / nu_new;
            post = NormalInvChisq(k_new, m_new, nu_new, s2_new);
        end

        % Log probability density function
        function val = logpdf(obj, mu, sigma2)
            scaledInvChiSq = ScaledInvChisq(obj.nu, obj.s2);
            val = logpdf(scaledInvChiSq, sigma2) + log(normpdf(mu, obj.m, sqrt(sigma2 / obj.k)));
        end
    end
end
