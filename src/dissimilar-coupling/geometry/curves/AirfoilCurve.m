classdef AirfoilCurve < SmoothCurveClosed
    % Joukowski airfoil: the image z = zeta + 1/zeta of the circle
    % zeta(t) = mu + (1 + regularized)*(1 - mu)*exp(2*pi*i*t), t in
    % [0, 1], with mu = mu_real + i*mu_imag.
    %
    %   curve = AirfoilCurve(mu_real, mu_imag, regularized)
    %
    % mu_real, mu_imag (defaults -1/10, 1/10) set the center of the
    % preimage circle; regularized (default 0) rounds the trailing edge:
    % 0 gives the exact airfoil with a cusp at t = 0, where the tangent
    % vanishes. mu (complex) and regularized are kept as properties.
    %
    % See also: SmoothCurveClosed
    properties
        mu
        regularized
    end

    methods
        function obj = AirfoilCurve(mu_real, mu_imag, regularized)
            arguments
                mu_real = -1/10
                mu_imag = 1/10
                % trailing-edge rounding: the preimage circle radius is
                % (1 + regularized)*|1 - mu|, enclosing the Joukowski
                % critical point zeta = 1; 0 = exact airfoil with a cusp
                % at t = 0, where the tangent vanishes
                regularized = 0
            end
            % zeta(t) = mu + (1 + regularized)*(1 - mu)*e^(2*pi*i*t),
            % z = zeta + 1/zeta:
            % with zeta = p + i q, p = mu_real + a cos(2 pi t) + b sin(2 pi t),
            % q = mu_imag + a sin(2 pi t) - b cos(2 pi t).
            a = (1 + regularized)*(1 - mu_real);
            b = (1 + regularized)*mu_imag;
            p = sprintf("%.17g + %.17g*cos(2*pi*t) + %.17g*sin(2*pi*t)", mu_real, a, b);
            q = sprintf("%.17g + %.17g*sin(2*pi*t) - %.17g*cos(2*pi*t)", mu_imag, a, b);

            x_t = sprintf("(%s) * (1 + 1/((%s)^2 + (%s)^2))", p, p, q);
            y_t = sprintf("(%s) * (1 - 1/((%s)^2 + (%s)^2))", q, p, q);

            obj = obj@SmoothCurveClosed(x_t, y_t, [0 1]);
            obj.mu = mu_real + 1i*mu_imag;
            obj.regularized = regularized;
        end
    end
end
