function [fluid_problem, solid_problem, problem_data, reference, sources] = lagrangian_march_problems( ...
    k, geometry, time, time_step, amplitude, options)
%LAGRANGIAN_MARCH_PROBLEMS  Manufactured transient FSI with a Lagrangian solid interface.
%
%   [fp, sp, pd, reference, sources] = lagrangian_march_problems(k, geometry, t, dt, A)
%
%   The solid is Lagrangian on its reference domain, below y = 0.5
%   (geometry = "flat") or inside the circle of radius 0.4 about
%   (0.5, 0.5) ("circle"); reference is that Sigma_0. Its exact
%   displacement u_S(X, t) carries the factor t^4 (start from rest), so
%   the exact interface is Sigma(t) = X + u_S(X, t). The fluid fields are
%   Eulerian in x, also with the factor t^4. Problems and fields at time t
%   come from FluidStructureCoupling.load_problems (densities 1, mu = 1,
%   E = 10, nu = 0.3 by default; pd.u_S, pd.u_t_S, pd.sigma_S, pd.u_F,
%   pd.sigma_F, ...).
%
%   solid_model = "hyperelastic" makes the solid a HyperelasticProblem
%   (material "svk", Saint Venant-Kirchhoff, or "neo_hookean"; its
%   default stabilization data) with the same exact u_S: pd.sigma_S is
%   then the first Piola stress P(I + grad u_S) and pd.f_S = rho_s u_S_tt
%   - div P (HyperelasticProblem.get_manufactured_solution), so that
%   sigma_S n_S below is the reference traction P N on Sigma_0.
%   solid_viscosity = mu_v (hyperelastic only) adds the Kelvin-Voigt
%   stress F S_v to the solid and to pd.sigma_S
%   (HyperelasticProblem, Data.solid_viscosity).
%
%   convection = true adds the Navier-Stokes convective term to the
%   fluid source, rho_f (u_F.grad) u_F at time t, for the Oseen fluid of
%   lagrangian_carved_march (fluid_model = "navier_stokes"); pd.f_F is
%   left as the Stokes source.
%
%   The interface data are built on the DISCRETE curve the march uses,
%   sources(curve) = {T1, T2} for a DeformedCurve over reference, with the
%   model of InterfaceProblem.fsi:
%     T1(x(s)) = J(s) sigma_F(x(s)) n_F(s) + sigma_S(X(s)) n_S(s),
%       J = |x'|/|X'|, n_F = -curve.unit_normal, n_S = reference.unit_normal,
%       evaluated at fluid points, s recovered by Gauss-Newton on the curve;
%     T2(X(s)) = u_F(x(s)) - d/dt u_S(X(s)), evaluated at solid points.
%   The exact fields then satisfy the coupling on that curve; the march's
%   geometry error is measured separately against pd.u_S on Sigma_0.
%
%   exact_geometry = true builds the same data on the EXACT curve
%   x(s) = X(s) + u_S(X(s), t) instead, whatever the curve argument: T1
%   takes J, n_F and the fluid stress point x(s) from it (s still
%   recovered from the fluid point on the discrete curve), T2 evaluates
%   u_F at x(s). The tangent is X'(s) + grad u_S(X(s)) X'(s), exact. The
%   data are then those of the continuous problem, and the errors of the
%   march include its geometry error (predictor, lag, reconstruction).
%
%See also: lagrangian_carved_march, solid_interface_displacement, DeformedCurve
arguments
    k (1,1) double
    geometry (1,1) string {mustBeMember(geometry, ["flat", "circle"])}
    time (1,1) double
    time_step (1,1) double
    amplitude (1,1) double
    options.uu_fluid (1,1) string = "t^4*cos(x + 2*y)*exp(x)/2"
    options.p_fluid (1,1) string = "t^4*sin(pi*x)*exp(y - 0.5)"
    options.rho_f = 1;  options.mu = 1;  options.rho_s = 1;  options.E = 10;  options.nu = 0.3;
    % flat only: u_S = 0 on the bottom y = 0 (the solid is clamped to the
    % floor and the domain stays the unit square; the factor 3y keeps the
    % interface displacement of the unclamped (1 + y) at y = 0.5); false
    % gives the unclamped data, whose bottom edge moves with the data
    options.clamped (1,1) logical = true
    options.solid_model (1,1) string {mustBeMember(options.solid_model, ["linear", "hyperelastic"])} = "linear"
    options.material (1,1) string {mustBeMember(options.material, ["svk", "neo_hookean"])} = "svk"
    options.solid_viscosity (1,1) double = 0
    % multiplies the fluid fields (uu_fluid, p_fluid); amplitude scales
    % the solid only, so fluid_scale = amplitude/A0 makes the whole
    % solution proportional to the amplitude
    options.fluid_scale (1,1) double = 1
    options.convection (1,1) logical = false
    options.exact_geometry (1,1) logical = false
end
if options.fluid_scale ~= 1
    scale = sprintf("%.17g*", options.fluid_scale);
    options.uu_fluid = scale + "(" + options.uu_fluid + ")";
    options.p_fluid = scale + "(" + options.p_fluid + ")";
end
a = sprintf("%.17g*t^4", amplitude);
if geometry == "flat"
    reference = FlatCurveHorizontal(0.5);
    profile = "(1 + y)";
    if options.clamped, profile = "(3*y)"; end
    u_solid = [a + "*sin(pi*x)*" + profile; a + "*cos(pi*x)*" + profile];
    parameter_of = @(P) P(1, :);
else
    reference = CircularCurve([0.5; 0.5], 0.4);
    u_solid = [a + "*(cos(pi*y) + x)"; a + "*(sin(pi*x) + y)"];
    parameter_of = @(P) mod(atan2(P(2, :) - 0.5, P(1, :) - 0.5)/(2*pi), 1);
end
[fluid_problem, solid_problem, ~, problem_data] = FluidStructureCoupling.load_problems( ...
    k, reference, time, time_step, options.uu_fluid, options.p_fluid, u_solid, ...
    options.rho_f, options.mu, options.rho_s, options.E, options.nu, 0, ...
    pressure_mean = "pinned");
if options.solid_model == "hyperelastic"
    material = options.material;
    if material == "svk", material = "saint_venant_kirchhoff"; end
    exact = HyperelasticProblem.get_manufactured_solution(u_solid, options.E, options.nu, ...
        0, options.rho_s, time, material, options.solid_viscosity);
    problem_data.sigma_S = exact.sigma;
    problem_data.f_S = exact.f;
    linear = solid_problem;
    solid_problem = HyperelasticProblem(k, material = material);
    for name = string(fieldnames(linear.Data)).'
        solid_problem.Data.(name) = linear.Data.(name);
    end
    solid_problem.Data.source_term = exact.f;
    solid_problem.Data.solid_viscosity = options.solid_viscosity;
end
if options.convection
    x = sym("x");  y = sym("y");
    stream = subs(str2sym(options.uu_fluid), sym("t"), time);
    u = [-diff(stream, y); diff(stream, x)];
    convective = options.rho_f*jacobian(u, [x, y])*u;
    f = fluid_problem.Data.source_term;
    for i = 1:2
        c = matlabFunction(convective(i), Vars = [x, y]);
        fluid_problem.Data.source_term{i} = @(x, y) f{i}(x, y) + c(x, y) + 0*x;
    end
end
pd = problem_data;
exact = [];
if options.exact_geometry
    exact = exact_curve(reference, u_solid, time, pd.u_S);
end
sources = @(curve) interface_sources(curve, reference, pd, parameter_of, exact);
end

function curve = exact_curve(reference, u_solid, time, u_S)
% x(s) = X(s) + u_S(X(s), t) with the analytic tangent (I + grad u_S) X'
x = sym("x");  y = sym("y");
grad_u = cell(2, 2);
for i = 1:2
    u = subs(str2sym(u_solid(i)), sym("t"), time);
    grad_u{i, 1} = matlabFunction(diff(u, x), Vars = [x, y]);
    grad_u{i, 2} = matlabFunction(diff(u, y), Vars = [x, y]);
end
d = @(s) on_points(u_S, reference.get_2d_points(s));
d_prime = @(s) gradient_times(grad_u, reference.get_2d_points(s), reference.tangent(s) + 0*s);
curve = DeformedCurve(reference, d, parametric = true, derivative = d_prime);
end

function v = on_points(u, P)
v = [u{1}(P(1, :), P(2, :)) + 0*P(1, :); u{2}(P(1, :), P(2, :)) + 0*P(1, :)];
end

function v = gradient_times(G, P, T)
% (grad u)(P) T, columnwise
v = zeros(2, size(P, 2));
for i = 1:2
    v(i, :) = (G{i, 1}(P(1, :), P(2, :)) + 0*P(1, :)) .* T(1, :) ...
        + (G{i, 2}(P(1, :), P(2, :)) + 0*P(1, :)) .* T(2, :);
end
end

function src = interface_sources(curve, reference, pd, parameter_of, exact)
% exact: [] (data on curve) or the exact curve sharing s with curve
on = curve;
if ~isempty(exact), on = exact; end
src = {{@(x, y) traction(x, y, 1); @(x, y) traction(x, y, 2)}, ...
       {@(x, y) kinematic(x, y, 1); @(x, y) kinematic(x, y, 2)}};

    function T = traction(x, y, i)
        sz = size(x);
        P = [x(:).'; y(:).'];
        t = parameter_of(P);
        for it = 1:20                            % Gauss-Newton for x(t) = P
            g = curve.tangent(t);
            step = sum((curve.get_2d_points(t) - P) .* g, 1) ./ sum(g.^2, 1);
            t = t - step;
            if max(abs(step)) < 1e-14, break; end
        end
        X = reference.get_2d_points(t);
        J = vecnorm(on.tangent(t), 2, 1) ./ vecnorm(reference.tangent(t) + 0*t, 2, 1);
        nF = -reshape(on.unit_normal(t), 2, []);
        nS = reshape(reference.unit_normal(t), 2, []) + 0*t;
        x = P;
        if ~isempty(exact), x = exact.get_2d_points(t); end
        T = J .* (pd.sigma_F{i, 1}(x(1, :), x(2, :)) .* nF(1, :) ...
                + pd.sigma_F{i, 2}(x(1, :), x(2, :)) .* nF(2, :)) ...
            + pd.sigma_S{i, 1}(X(1, :), X(2, :)) .* nS(1, :) ...
            + pd.sigma_S{i, 2}(X(1, :), X(2, :)) .* nS(2, :);
        T = reshape(T, sz);
    end

    function T = kinematic(x, y, i)
        sz = size(x);
        X = [x(:).'; y(:).'];
        xs = on.get_2d_points(parameter_of(X));
        T = pd.u_F{i}(xs(1, :), xs(2, :)) - pd.u_t_S{i}(X(1, :), X(2, :));
        T = reshape(T, sz);
    end
end
