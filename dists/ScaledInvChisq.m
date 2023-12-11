classdef ScaledInvChisq
    properties
        nu    % Degrees of freedom
        s2    % Scale parameter
    end

    methods
        % Constructor
        function obj = ScaledInvChisq(nu, s2)
            obj.nu = nu;
            obj.s2 = s2;
        end

        % Log probability density function
        function val = logpdf(obj, theta)
            hnu = 0.5 * obj.nu;
            val = hnu * log(hnu) - gammaln(hnu) + hnu * log(obj.s2) - (hnu + 1) * log(theta) - hnu * obj.s2 / theta;
        end

        % Random number generation
        function r = rand(obj)
            r = obj.nu * obj.s2 / chi2rnd(obj.nu);
        end
    end
end
