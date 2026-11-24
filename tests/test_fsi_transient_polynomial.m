% Transient FSI time stepping on a solution both integrators reproduce
% exactly, on a flat Sigma with fitted meshes at gap 0 and gap 0.02:
%
%   solid    u_S = A(x,y) * g(t),   g(t) = 1 + t + t^2/2,   g'' = 1
%   fluid    u_F = B(x,y) * (1 + t)
%
% Newmark(beta = 1/4, gamma = 1/2) is exact for constant acceleration and
% Crank-Nicolson for velocity linear in t; A and B are degree <= k and the
% flat interface is represented exactly. The march starts from
% FluidStructureCoupling.initial_state (which reproduces the projections
% of the exact accelerations for these data) and carries the state
% returned by FluidStructureCoupling.solve; the T2 history is transported
% onto Sigma, so both gaps are asserted: the error must stay below 1e-10
% at every step.
%
%See also: test_trace_extension, FluidStructureCoupling.solve,
%          FluidStructureCoupling.initial_state
clc; clear;

k    = 2;
h    = 1/8;
dt   = 0.1;
n_steps = 4;
E    = 1.0;
nu   = 0.25;
visc = 0.4;
rho_f = 1.0;      % both densities must be > 0 or the problems are steady
rho_s = 1.0;

tol = 1e-10;

% Degree <= k in space, with the time factors above.
uu_fluid = "x*y*(1 + t)";            % u_F = curl of this: linear in x and t
p_fluid  = "(x + y)*(1 + t)";
u_solid  = ["x*(1 + t + t^2/2)"; "y*(1 + t + t^2/2)"];

gaps = [0, 0.02];
interface = FlatCurveHorizontal();

n_fail = 0;
for gi = 1:numel(gaps)
    g = gaps(gi);

    [fluid_mesh, solid_mesh] = fitted_split(interface, h, g);
    fprintf("\ngap = %g\n", g);

    % t = 0 data, and the initial state
    [pde_fluid, pde_solid, sources, d0] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
        uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
    pde_solid.Data.time_scheme = "newmark";
    % Crank-Nicolson from the first step: with the consistent a_f^0 of the
    % initial state the trapezoid rule is already exact, a Backward-Euler
    % startup step would not be.
    pde_fluid.Data.time_scheme = "cranknicolson";
    state = FluidStructureCoupling.initial_state(pde_fluid, pde_solid, fluid_mesh, solid_mesh, ...
        interface, sources{:}, source_dirichlet_rate = d0.T_2_t);

    for n = 1:n_steps
        t = n*dt;
        [pde_fluid, pde_solid, sources, pd] = ...
            FluidStructureCoupling.load_problems(k, interface, t, dt, ...
                uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
        pde_solid.Data.time_scheme = "newmark";
        pde_fluid.Data.time_scheme = "cranknicolson";

        [sol_fluid, sol_solid, state] = FluidStructureCoupling.solve( ...
            pde_fluid, pde_solid, fluid_mesh, solid_mesh, interface, ...
            sources{:}, "newmark", quiet = true, state = state);

        err_s = sol_solid.compute_l2_error( ...
            {pd.sigma_S, pd.u_S, pd.gamma_S}, {pd.u_S});
        err_f = sol_fluid.compute_l2_error( ...
            {pd.sigma_F, pd.u_F, pd.gamma_F}, {pd.u_F});
        worst = max([err_s(:); err_f(:)]);

        if worst > tol
            fprintf("  step %d (t = %.2f): NOT exact, worst L2 = %.3e  " + ...
                "[solid %.2e, fluid %.2e]\n", n, t, worst, ...
                max(err_s(:)), max(err_f(:)));
            n_fail = n_fail + 1;
        else
            fprintf("  step %d (t = %.2f): ok (worst L2 = %.2e)\n", n, t, worst);
        end
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_transient_polynomial:failed", "%d check(s) failed", n_fail);
end
