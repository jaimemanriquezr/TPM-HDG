% T2 history transport: the coupling transports the T2 history onto Sigma
% itself (InterfaceProblem.solve, from state.solid), which
% makes the gapped transient FSI coupling exact on a polynomial solution
% (flat Sigma, fitted meshes retreated by 0.02, one Newmark /
% Crank-Nicolson step from FluidStructureCoupling.initial_state at several
% dt).
%
% The history is analytic: with beta = 1/4, gamma = 1/2 and g(0) = g'(0) = 1
% it is A*(1 + dt/2), whose gradient is (1 + dt/2)*I for A = [x; y]. Its
% extension, built by the same Legendre projection trace_extension_coefficients
% uses, must equal the one the coupling forms (recomputed here from the
% state's (sigma, gamma) triple, and from the (sigma, gamma) the coupled
% problem stores) to
% 1e-12, and the step must be exact (1e-10). Transport is not optional, so
% there is no arm without it.
%
%See also: test_fsi_transient_polynomial, test_trace_extension
clc; clear;

k = 2; h = 1/8; E = 1.0; nu = 0.25; visc = 0.4; rho_f = 1.0; rho_s = 1.0;
beta = 1/4; gamma = 1/2;
gap = 0.02;
dts = [0.2, 0.1, 0.05];

tol = 1e-10;
tol_extension = 1e-12;

uu_fluid = "x*y*(1 + t)";
p_fluid  = "(x + y)*(1 + t)";
u_solid  = ["x*(1 + t + t^2/2)"; "y*(1 + t + t^2/2)"];

interface = FlatCurveHorizontal();
[fluid_mesh, solid_mesh] = fitted_split(interface, h, gap);
[~, part_S] = interface.get_partitions(fluid_mesh, solid_mesh, ...
    quadrature_patch = k + 1, quadrature_path = k + 2);

d1      = k + 1;
s_quad  = size(part_S.quadrature.local, 2);
[~, wS] = lgwt(s_quad, 0, 1);
n_patch = size(part_S.faces, 3);

n_fail = 0;
for dt = dts
    [pde_f0, pde_s0, src0, d0] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
        uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
    pde_s0.Data.time_scheme = "newmark";
    pde_f0.Data.time_scheme = "cranknicolson";
    state = FluidStructureCoupling.initial_state(pde_f0, pde_s0, fluid_mesh, solid_mesh, ...
        interface, src0{:}, source_dirichlet_rate = d0.T_2_t);

    [pde_f, pde_s, src, pd] = FluidStructureCoupling.load_problems(k, interface, dt, dt, ...
        uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
    pde_s.Data.time_scheme = "newmark";
    pde_f.Data.time_scheme = "cranknicolson";

    % the analytic extension
    grad_h0 = (1 + dt/2)*eye(2);
    hx = zeros(2*d1, n_patch);
    for e = 1:n_patch
        m = part_S.quadrature.interface(:, :, e) - part_S.quadrature.physical(:, :, e);
        l = legendre_basis(k, part_S.quadrature.local(1, :, e));
        c = (l' * (wS .* l)) \ (l' * (wS .* (grad_h0 * m).'));
        hx(:, e) = [c(:, 1); c(:, 2)];
    end
    % the coupling's, recomputed from the state and read back
    s = state.solid;
    c_v = dt*(1 - beta/gamma);
    c_a = dt^2*((1/2 - beta) - beta*(1 - gamma)/gamma);
    hx_state = trace_extension_coefficients(pde_s, solid_mesh, part_S, ...
        s.stress + c_v*s.stress_velocity + c_a*s.stress_acceleration, ...
        s.rotation + c_v*s.rotation_velocity + c_a*s.rotation_acceleration, Order = k);
    coupled = FluidStructureCoupling.interface_problem(pde_f, pde_s, fluid_mesh, solid_mesh, ...
        interface, src{:}, "newmark", state = state);
    d_ext = max(abs(hx_state - hx), [], "all");
    same = isequal(hx_state, trace_extension_coefficients(pde_s, solid_mesh, part_S, ...
        coupled.coupling.history_stress, coupled.coupling.history_rotation, Order = k));

    [sol_f, sol_s] = FluidStructureCoupling.solve(pde_f, pde_s, ...
        fluid_mesh, solid_mesh, interface, src{:}, "newmark", quiet = true, state = state);

    es = sol_s.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S}, {pd.u_S});
    ef = sol_f.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F}, {pd.u_F});
    worst = max([es(:); ef(:)]);

    fprintf("dt=%-6g |Ext(state) - Ext(analytic)| %.2e  coupling's == recomputed %d   " + ...
        "solid %.3e   fluid %.3e\n", dt, d_ext, same, max(es(:)), max(ef(:)));

    if d_ext > tol_extension || ~same
        fprintf("  FAILED: the coupling's extension is not the analytic one\n");
        n_fail = n_fail + 1;
    end
    if worst > tol
        fprintf("  FAILED: the transported step is not exact\n");
        n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_transient_transport:failed", "%d check(s) failed", n_fail);
end
