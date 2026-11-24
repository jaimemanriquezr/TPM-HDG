% Gapped transient FSI on carved meshes, with an exact reference: a solution
% both integrators reproduce exactly (constant solid acceleration for
% Newmark(1/4,1/2), linear fluid velocity for Crank-Nicolson, degree <= k in
% space), carved out of a background grid, on a flat and a curved Sigma.
% The march starts from FluidStructureCoupling.initial_state and carries
% the state FluidStructureCoupling.solve returns; the coupling transports
% the T2 history onto Sigma itself (there is no arm without transport).
%
% Only the flat arm is asserted: the error stays at round-off (tol 1e-8).
% The curved arm is reported, not failed.
%
% Keep N odd so y = 0.5 is not a grid line (else the carve is fitted, gap 0),
% and keep the "left" diagonal ("right" makes the coupled solve near-singular).
%
%See also: test_fsi_transient_transport, test_fsi_transient_polynomial
clc; clear;

k = 2; E = 1.0; nu = 0.25; visc = 0.4; rho_f = 1.0; rho_s = 1.0;
N = 9; dt = 1e-3; n_steps = 12;

tol_exact = 1e-8;    % flat: exact up to the recursion's 1/(beta dt^2)

uu_fluid = "x*y*(1 + t)";
p_fluid  = "(x + y)*(1 + t)";
u_solid  = ["x*(1 + t + t^2/2)"; "y*(1 + t + t^2/2)"];

n_fail = 0;
for is_curved = [false, true]
    if is_curved
        interface = SmoothCurveHorizontal();
    else
        interface = FlatCurveHorizontal();
    end
    [fluid_mesh, solid_mesh] = StructuredSimplicialGrid(N, N, "left") ...
        .split(interface, cut_elements = "discard");

    [pde_f, pde_s, src, d0] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
        uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
    pde_s.Data.time_scheme = "newmark";
    pde_f.Data.time_scheme = "cranknicolson";
    state = FluidStructureCoupling.initial_state(pde_f, pde_s, fluid_mesh, solid_mesh, ...
        interface, src{:}, source_dirichlet_rate = d0.T_2_t);
    worst = NaN; blew = 0;

    for n = 1:n_steps
        t = n*dt;
        [pde_f, pde_s, src, pd] = FluidStructureCoupling.load_problems(k, interface, t, dt, ...
            uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
        pde_s.Data.time_scheme = "newmark";
        pde_f.Data.time_scheme = "cranknicolson";

        [sol_f, sol_s, state] = FluidStructureCoupling.solve(pde_f, pde_s, fluid_mesh, solid_mesh, ...
            interface, src{:}, "newmark", quiet = true, state = state);

        es = sol_s.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S}, {pd.u_S});
        ef = sol_f.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F}, {pd.u_F});
        worst = max([es(:); ef(:)]);
        if ~isfinite(worst) || worst > 1e6, blew = n; break; end
    end

    if blew
        fprintf("curved=%d   DIVERGED at step %d\n", is_curved, blew);
    else
        fprintf("curved=%d   final %.3e\n", is_curved, worst);
    end

    if ~is_curved && (blew || worst > tol_exact)
        % the flat arm is the assertion
        fprintf("  FAILED: flat carve is not exact\n");
        n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_curved_carve:failed", "%d check(s) failed", n_fail);
end
