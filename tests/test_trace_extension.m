% Consistency of the TPM trace extension on polynomial displacement fields.
% trace_extension_coefficients returns
%
%     Ext(u) := int_0^1 grad u(x + xi*m) . m  dxi,     m = xbar - x,
%
% the term the T2 kinematic row adds to a raw skeleton trace to land it on
% Sigma, with grad u reconstructed from the PEERS variables as
% c1*sigma - c2*tr(sigma)*I + gamma (an identity for the exact fields). For
% affine grad u the path integral has the closed form grad u(x + m/2) . m.
%
% Per case, from an actual coupled solve: check 1, the solid polynomial is
% reproduced to 1e-10; check 2, Ext matches the closed form to 1e-10. Both
% flat runs are asserted (gap = 0, the inert control, and gap = 0.02, where
% m ~= 0); curved runs are reported.
%
%See also: trace_extension_coefficients, FluidStructureCoupling.solve
clc; clear;

k    = 2;
h    = 1/8;
gap  = 0.02;      % retreat both polygons from Sigma so that m ~= 0
E    = 1.0;
nu   = 0.25;
visc = 0.4;

tol_poly = 1e-10;  % exact reproduction of a polynomial
tol_ext  = 1e-10;  % agreement with the closed-form extension

% The fluid fields must also be of degree <= k: the sides are coupled, so a
% fluid solution the space cannot represent leaks its discretization error
% into the solid.
uu_fluid = "x*y";
p_fluid  = "x + y";

% gap = 0 is the control: both meshes sit on Sigma, m = 0, and the transfer
% is inert.
gaps = [0, gap];

% fitted_split meshes each side as a polygon, so a curved Sigma carries an
% O(h^2) geometry error whatever the field; only a flat Sigma can be held
% to exactness.
interfaces = {FlatCurveHorizontal(), SmoothCurveHorizontal()};

% Each case: solid displacement (symbolic strings) and its exact gradient as
% a function of position, grad(x,y) -> 2x2 with rows [du_i/dx, du_i/dy].
cases = { ...
    "rigid translation", ["1";        "2"],  @(x, y) zeros(2, 2); ...
    "uniform stretch",   ["x";        "y"],  @(x, y) eye(2); ...
    "pure rotation",     ["y";       "-x"],  @(x, y) [0 1; -1 0]; ...
    "quadratic",         ["x^2";      "0"],  @(x, y) [2*x 0; 0 0]};

n_fail = 0;
for ii = 1:numel(interfaces)
interface = interfaces{ii};
is_flat   = isa(interface, "FlatCurveHorizontal");
for gi = 1:numel(gaps)
g = gaps(gi);

% Strict where the geometry is exact (flat Sigma); curved configurations
% are reported, not failed.
strict = is_flat;

[upper_mesh, lower_mesh] = fitted_split(interface, h, g);

% mesh_2 / lower is the solid side, the one whose history the T2 row carries.
[~, partition] = interface.get_partitions(upper_mesh, lower_mesh, ...
    quadrature_patch = k + 1, quadrature_path = k + 2);

if strict
    fprintf("\n%s, gap = %g   [strict]\n", class(interface), g);
else
    fprintf("\n%s, gap = %g\n", class(interface), g);
end

for c = 1:size(cases, 1)
    name    = cases{c, 1};
    u_solid = cases{c, 2};
    grad_fn = cases{c, 3};

    [fluid_problem, solid_problem, interface_sources, problem_data] = ...
        FluidStructureCoupling.load_problems(k, interface, 0, 1e-3, ...
            uu_fluid, p_fluid, u_solid, ...
            0, visc, 0, E, nu, 0);

    [~, solid_solution] = FluidStructureCoupling.solve( ...
        fluid_problem, solid_problem, upper_mesh, lower_mesh, interface, ...
        interface_sources{:}, quiet = true);

    % Check 1: the gapped TPM solve reproduces the polynomial exactly.
    err = solid_solution.compute_l2_error( ...
        {problem_data.sigma_S, problem_data.u_S, problem_data.gamma_S}, ...
        {problem_data.u_S});
    err_u = max(err(:));
    if err_u > tol_poly
        fprintf("  [%s] polynomial NOT reproduced: L2 = %.3e\n", name, err_u);
        n_fail = n_fail + strict;
    end

    % Check 2: Ext against the closed form.
    sigma_dof = solid_solution.Unknowns.field(1).dof;
    gamma_dof = solid_solution.Unknowns.field(3).dof;
    if iscell(gamma_dof), gamma_dof = cat(1, gamma_dof{:}); end

    ext = trace_extension_coefficients(solid_problem, lower_mesh, ...
        partition, sigma_dof, gamma_dof, Order = k);

    d1      = k + 1;
    s_quad  = size(partition.quadrature.local, 2);
    [~, wS] = lgwt(s_quad, 0, 1);
    n_patch = size(partition.faces, 3);

    worst = 0;
    for e = 1:n_patch
        xq = partition.quadrature.physical(:, :, e);              % 2 x s
        m  = partition.quadrature.interface(:, :, e) - xq;        % 2 x s

        % int_0^1 grad u(x + xi m) . m dxi = grad u(x + m/2) . m, exact for
        % affine grad u (covers every case above).
        ext_exact = zeros(s_quad, 2);
        for q = 1:s_quad
            xm = xq(:, q) + 0.5*m(:, q);
            ext_exact(q, :) = (grad_fn(xm(1), xm(2)) * m(:, q)).';
        end

        % Same L2 projection onto the face Legendre basis the helper uses, so
        % only the reconstruction and the path integral are being compared.
        l = legendre_basis(k, partition.quadrature.local(1, :, e));
        coeffs = (l' * (wS .* l)) \ (l' * (wS .* ext_exact));
        worst = max(worst, max(abs(ext(:, e) - [coeffs(:, 1); coeffs(:, 2)])));
    end

    if worst > tol_ext
        fprintf("  [%s] Ext MISMATCH: max |code - exact| = %.3e\n", name, worst);
        n_fail = n_fail + strict;
    else
        fprintf("  [%s] ok (L2 = %.2e, Ext err = %.2e)\n", name, err_u, worst);
    end
end
end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_trace_extension:failed", "%d check(s) failed", n_fail);
end
