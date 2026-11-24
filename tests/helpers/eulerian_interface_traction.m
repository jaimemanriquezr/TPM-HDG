function T1 = eulerian_interface_traction(curve, sigma_F, sigma_S, options)
%EULERIAN_INTERFACE_TRACTION  Traction datum of the Eulerian FSI coupling on a discrete curve.
%
%   T1 = eulerian_interface_traction(curve, sigma_F, sigma_S, n_samples = 2000)
%
%   The datum T1 = sigma_F n_F + sigma_S n_S of the traction balance
%   (load_problems' T_1) with both Cauchy stresses at the same point x and
%   the normals of curve: n_S = n, n_F = -n, n = curve.unit_normal
%   (upward over a graph, outward around a closed curve, i.e. from the
%   solid into the fluid),
%       T1{i}(x, y) = sum_j (sigma_S{i,j}(x) - sigma_F{i,j}(x)) n_j(t(x)),
%   t(x) the parameter of the closest point of curve to x: the nearest of
%   n_samples points of curve evenly spaced in its parameter, polished by
%   Gauss-Newton steps on |x(t) - x|^2 (at most 20, to a step below
%   1e-14). The exact fields then satisfy the traction balance on curve,
%   whatever curve is; on the exact Sigma(t) T1 equals load_problems' T_1.
%   sigma_F, sigma_S: 2 x 2 cells of handles (x, y) (pd.sigma_F,
%   pd.sigma_S). The parameters of the last points are reused when the
%   same points are queried again (the two components).
%
%See also: eulerian_solid_march, translating_flat_problems
arguments
    curve
    sigma_F cell
    sigma_S cell
    options.n_samples (1,1) double {mustBeInteger, mustBePositive} = 2000
end
range = curve.t_range;
periodic = is_periodic(curve);
if periodic
    t_grid = range(1) + diff(range)*(0:options.n_samples - 1)/options.n_samples;
else
    t_grid = linspace(range(1), range(2), options.n_samples);
end
P_grid = reshape(curve.get_2d_points(t_grid), 2, []);
P_last = [];
t_last = [];
T1 = {@(x, y) traction(x, y, 1); @(x, y) traction(x, y, 2)};

    function T = traction(x, y, i)
        sz = size(x);
        P = [x(:).'; y(:).'];
        t = parameter(P);
        n = reshape(curve.unit_normal(t), 2, []) + 0*t;
        T = (sigma_S{i, 1}(P(1, :), P(2, :)) - sigma_F{i, 1}(P(1, :), P(2, :))).*n(1, :) ...
            + (sigma_S{i, 2}(P(1, :), P(2, :)) - sigma_F{i, 2}(P(1, :), P(2, :))).*n(2, :);
        T = reshape(T, sz);
    end

    function t = parameter(P)
        if isequal(P, P_last)
            t = t_last;
            return
        end
        t = zeros(1, size(P, 2));
        for j = 1:size(P, 2)
            [~, nearest] = min(sum((P_grid - P(:, j)).^2, 1));
            t(j) = t_grid(nearest);
        end
        for it = 1:20
            g = reshape(curve.tangent(t), 2, []) + 0*t;
            step = sum((reshape(curve.get_2d_points(t), 2, []) - P) .* g, 1) ./ sum(g.^2, 1);
            t = t - step;
            if periodic
                t = range(1) + mod(t - range(1), diff(range));
            else
                t = min(max(t, range(1)), range(2));
            end
            if max(abs(step)) < 1e-14, break; end
        end
        [P_last, t_last] = deal(P, t);
    end
end
