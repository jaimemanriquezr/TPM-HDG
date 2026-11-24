function [fluid_problem, solid_problem, sources, problem_data, curves, mesher] = ...
    lagrangian_pair_problems(k, geometry, options)
%LAGRANGIAN_PAIR_PROBLEMS  Steady manufactured FSI on a per-side curve pair.
%
%   [fp, sp, sources, pd, curves, mesher] = lagrangian_pair_problems(k, geometry)
%
%   The solid lives on its reference domain, below the line y = 0.5
%   (geometry = "flat") or inside the circle of radius 0.4 about
%   (0.5, 0.5) ("circle"): Sigma_0, s -> X(s), the reference curve's own
%   parametrization. The fluid lies above / outside the deformed curve
%   Sigma = (I + d)(Sigma_0), x(s) = X(s) + d(X(s)), with the prescribed
%   smooth displacement
%     flat    d(X) = A [sin(pi X1); cos(pi X1)],        A = 0.05
%     circle  d(X) = A [cos(3 pi X2); sin(2 pi X1)],    A = 0.03
%   curves = {DeformedCurve(Sigma_0, d), Sigma_0}, the fsi pair.
%
%   Exact fields on R^2 from FluidStructureCoupling.load_problems (steady,
%   time_step 1, pinned pressure). The interface data follow the model of
%   InterfaceProblem.fsi, with x'(s), the normals and the Nanson factor
%   J = |x'|/|X'| derived symbolically (independent of DeformedCurve's
%   finite differences):
%     T1(x(s)) = J(s) sigma_F(x(s)) n_F(s) + sigma_S(X(s)) n_S(s),
%       n_F the fluid's outward normal of Sigma, n_S the solid's of
%       Sigma_0; evaluated at fluid points, s recovered by Gauss-Newton;
%     T2(X(s)) = u_F(x(s)) - u_S(X(s))/dt, evaluated at solid points.
%
%   mesher(h) returns [fluid_mesh, solid_mesh]:
%     flat    solid: structured on [0,1] x [0,0.5], N = 1/h cells across;
%             fluid: structured 3N/2 cells across, mapped by
%             (xi, eta) -> (1 - eta) x(xi) + eta (xi, 1), so its
%             interface vertices lie on Sigma;
%     circle  fitted_split of Sigma_0 (solid, inner, size h) and of Sigma
%             (fluid, outer, size 3h/4), interface vertices about one
%             element size apart, on the curves.
%
%See also: DeformedCurve, InterfaceProblem.fsi, FluidStructureCoupling.load_problems
arguments
    k (1,1) double
    geometry (1,1) string {mustBeMember(geometry, ["flat", "circle"])}
    options.uu_fluid (1,1) string = "cos(x + 2*y)*exp(x)/2"
    options.p_fluid (1,1) string = "sin(pi*x)*exp(y - 0.5)"
    options.u_solid (2,1) string = ["0.5*sin(2*x + y)"; "0.5*cos(x - y)"]
    options.amplitude double {mustBeScalarOrEmpty} = []
    % false drops J from T1 (a negative control: the data then miss the
    % model's Nanson factor)
    options.nanson (1,1) logical = true
end
dt = 1.0;
s = sym("s", "real");  a = sym("a", "real");  b = sym("b", "real");
if geometry == "flat"
    reference = FlatCurveHorizontal(0.5);
    X = [s; sym(0.5)];
    A = 0.05;  d = [sin(pi*a); cos(pi*a)];
    sense = 1;                                   % normal = tangent rotated +90
    guess = @(x, y) x;
else
    reference = CircularCurve([0.5; 0.5], 0.4);
    X = [0.5 + 0.4*cos(2*pi*s); 0.5 + 0.4*sin(2*pi*s)];
    A = 0.03;  d = [cos(3*pi*b); sin(2*pi*a)];
    sense = -1;                                  % outward, rotated -90
    guess = @(x, y) mod(atan2(y - 0.5, x - 0.5)/(2*pi), 1);
end
if ~isempty(options.amplitude), A = options.amplitude; end
d = A*d;
x_of_s = X + subs(d, [a, b], X.');
xp = diff(x_of_s, s);
Xp = diff(X, s);
rot = @(v) sense*[-v(2); v(1)]/sqrt(v(1)^2 + v(2)^2);
vec = @(e) vectorized(e, s);
x_fun   = vec(x_of_s);
xp_fun  = vec(xp);
X_fun   = vec(X);
J_fun   = vec(sqrt(xp.'*xp)/sqrt(Xp.'*Xp));
nF_fun  = vec(-rot(xp));                         % fluid outward normal of Sigma
nS_fun  = vec(rot(Xp));                          % solid outward normal of Sigma_0
d_fun   = arrayfun(@(c) matlabFunction(c, Vars = [a, b]), d, UniformOutput = false);

displacement = @(P) [d_fun{1}(P(1, :), P(2, :)) + 0*P(1, :); d_fun{2}(P(1, :), P(2, :)) + 0*P(1, :)];
curves = {DeformedCurve(reference, displacement), reference};

[fluid_problem, solid_problem, ~, problem_data] = FluidStructureCoupling.load_problems( ...
    k, reference, 0.0, dt, options.uu_fluid, options.p_fluid, options.u_solid, ...
    0.0, 1.0, 0.0, 1.0, 0.25, 0.0, pressure_mean = "pinned");
pd = problem_data;

    function T = traction(x, y, i)
        sz = size(x);
        P = [x(:).'; y(:).'];
        t = guess(P(1, :), P(2, :));
        for it = 1:30                            % Gauss-Newton for x(t) = P
            r = x_fun(t) - P;
            g = xp_fun(t);
            step = sum(r .* g, 1) ./ sum(g.^2, 1);
            t = t - step;
            if max(abs(step)) < 1e-15, break; end
        end
        assert(max(vecnorm(x_fun(t) - P)) < 1e-12, "lagrangian_pair_problems:inverse", ...
            "Gauss-Newton did not recover the parameter of a point of Sigma.");
        Xs = X_fun(t);  J = J_fun(t).^options.nanson;  nF = nF_fun(t);  nS = nS_fun(t);
        T = J .* (pd.sigma_F{i, 1}(P(1, :), P(2, :)) .* nF(1, :) ...
                + pd.sigma_F{i, 2}(P(1, :), P(2, :)) .* nF(2, :)) ...
            + pd.sigma_S{i, 1}(Xs(1, :), Xs(2, :)) .* nS(1, :) ...
            + pd.sigma_S{i, 2}(Xs(1, :), Xs(2, :)) .* nS(2, :);
        T = reshape(T, sz);
    end

    function T = kinematic(x, y, i)
        sz = size(x);
        P = [x(:).'; y(:).'];
        xs = P + displacement(P);
        T = pd.u_F{i}(xs(1, :), xs(2, :)) - pd.u_S{i}(P(1, :), P(2, :))/dt;
        T = reshape(T, sz);
    end

sources = {{@(x, y) traction(x, y, 1); @(x, y) traction(x, y, 2)}, ...
           {@(x, y) kinematic(x, y, 1); @(x, y) kinematic(x, y, 2)}};
problem_data.T_1 = sources{1};
problem_data.T_2 = sources{2};
problem_data.nanson = J_fun;
problem_data.tangent = xp_fun;
problem_data.normal_fluid = nF_fun;

if geometry == "flat"
    mesher = @(h) flat_meshes(h, x_fun);
else
    mesher = @(h) circle_meshes(h, curves);
end
end

function f = vectorized(e, s)
% a column of expressions in s as a handle on a row of parameters, 2 x n
g = arrayfun(@(c) matlabFunction(c, Vars = s), e, UniformOutput = false);
f = @(t) cell2mat(cellfun(@(gi) gi(t) + 0*t, g, UniformOutput = false));
end

function [fluid_mesh, solid_mesh] = flat_meshes(h, x_fun)
N = round(1/h);
solid_mesh = StructuredSimplicialGrid(N, N/2, 0, 1, 0, 0.5);
solid_mesh = apply_tag(solid_mesh, "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
NF = round(1.5*N);
ref = StructuredSimplicialGrid(NF, round(NF/2));
xi = ref.coordinates(1, :);  eta = ref.coordinates(2, :);
P = (1 - eta) .* x_fun(xi) + eta .* [xi; ones(size(xi))];
fluid_mesh = SimplicialGrid(P, ref.elements, mesh_generator = "array");
bottom = P(:, eta == 0);
fluid_mesh = apply_tag(fluid_mesh, "Interface", @(Q) ismember(Q.', bottom.', "rows").');
end

function [fluid_mesh, solid_mesh] = circle_meshes(h, curves)
% interface vertices about one element size apart (circumference 0.8 pi)
points = @(h) round(0.8*pi/h);
[~, solid_mesh] = fitted_split(curves{2}, h, interface_points = points);
fluid_mesh = fitted_split(curves{1}, 0.75*h, interface_points = points);
end
