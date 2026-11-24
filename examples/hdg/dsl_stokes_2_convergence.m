function dsl_stokes_2_convergence(k, max_iterations, stream_function, pressure, viscosity, stress_space)
% Stokes in the stress-velocity-vorticity form written in the form language
% and solved on the declarative path, with a manufactured solution on
% structured meshes of the unit square (zero-mean pressure; stress in PEERS
% or in plain P_k^(2x2), i.e. PEERS without its bubbles):
%   c1 sigma - c2 tr(sigma) I - grad u + gamma J = 0,  -div sigma = f,
%   (sigma, eta J) = 0,  u = g  on the boundary,
% c1 = 1/(2mu), c2 = 1/(4mu), J = [0 1; -1 0], flux stabilization mu*tau.
% Prints L2 errors of sigma_h, u_h, gamma_h, u_hat and EOC ~ k+1.
arguments
    k               = 2
    max_iterations  = 4
    stream_function = "x^3*y^3"
    pressure        = "x*y"
    viscosity       = 1e-3
    stress_space    (1,1) string {mustBeMember(stress_space, ["peers", "pk"])} = "peers"
end
mu = viscosity;
c1 = 1/(2*mu);  c2 = 1/(4*mu);
J  = [0, 1; -1, 0];
tau = 1;

data = IncompressibleElasticityProblem.get_manufactured_solution(stream_function, pressure, mu);
g = Datum(data.u);
f = Datum(data.f);

Pk = DubinerBasis(k);
if stress_space == "peers"
    S_k = PEERSBasis(k);
else
    S_k = Pk.to_tensor();
end
% local: sigma, u, gamma, xi (trace multiplier); global: u_hat, rho
% (per-element pressure mean), mu_c (zero-mean pressure multiplier)
spaces = [FESpace([S_k, Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
          FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])];

[error_sigma, error_u, error_gamma, error_uhat, h_list, N_list] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    N    = 2^i;
    mesh = StructuredSimplicialGrid(N, N).add_stabilization(tau);

    problem = FEProblem(Spaces = spaces, Order = k);
    problem.Mesh = mesh;
    problem.BoundaryConditions = WeakEssentialBoundaryConditions(g, "Dirichlet");

    [sg, u, gm, xi, u_hat, rho, mu_c] = get_trial_functions(problem.Space);
    [T, v, e, z, v_hat, s_rho, s_mu]  = get_test_functions(problem.Space);
    dx = Measure(preset = "dx");
    dS = Measure(preset = "dS");
    n_ = ExteriorNormal();
    tau_ = Stabilization();

    % local equations: constitutive, momentum, weak symmetry, trace = rho
    problem.Operator = AffineOperator( ...
          dot(sg, c1*T)*dx - trace(sg)*(c2*trace(T))*dx ...
        + dot(u, divergence(T))*dx + dot(gm*J, T)*dx + xi*trace(T)*dx ...
        - dot(u_hat*n_.', T)*dS ...
        - dot(divergence(sg), v)*dx + dot((mu*tau_)*u, v)*dS ...
        - dot((mu*tau_)*u_hat, v)*dS + dot(f, v)*dx ...
        + dot(sg, e*J)*dx ...
        + trace(sg)*z*dx - rho*z*dx);
    % global equations: normal-flux continuity, <u_hat.n, 1>_dK = 0,
    % zero-mean pressure
    problem.Operator = problem.Operator + ( ...
          dot(sg*n_, v_hat)*dS - dot((mu*tau_)*u, v_hat)*dS ...
        + dot((mu*tau_)*u_hat, v_hat)*dS ...
        + dot(u_hat, n_)*s_rho*dS ...
        + mu_c*s_rho*dx + rho*s_mu*dx);

    solution = solve(problem, mesh, k);
    % per field: the multiplier spaces (xi, rho, mu_c) have no exact solution
    fields = solution.Unknowns.field;

    N_list(i)       = N;
    h_list(i)       = sqrt(2) / N;
    error_sigma(i)  = compute_l2_error(fields(1), mesh, Datum(data.sigma));
    error_u(i)      = compute_l2_error(fields(2), mesh, g);
    error_gamma(i)  = compute_l2_error(fields(3), mesh, Datum(data.gamma));
    error_uhat(i)   = compute_l2_error(solution.Unknowns.hybrid(1), mesh, g);
end

r_sigma = [NaN; compute_eoc(error_sigma, h_list)];
r_u     = [NaN; compute_eoc(error_u,     h_list)];
r_gamma = [NaN; compute_eoc(error_gamma, h_list)];
r_uhat  = [NaN; compute_eoc(error_uhat,  h_list)];

hdr = "   N  |     h    |  e(σ_h)   | r(σ_h) |   e(u_h)  | r(u_h) |  e(γ_h)   | r(γ_h) |  e(û_h)   | r(û_h)";
fprintf("\n%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
for i = 1:max_iterations
    if i == 1
        fprintf(" %4d | %8.2e | %9.3e |   -    | %9.3e |   -    | %9.3e |   -    | %9.3e |   -\n", ...
            N_list(i), h_list(i), error_sigma(i), error_u(i), error_gamma(i), error_uhat(i));
    else
        fprintf(" %4d | %8.2e | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f\n", ...
            N_list(i), h_list(i), error_sigma(i), r_sigma(i), error_u(i), r_u(i), error_gamma(i), r_gamma(i), error_uhat(i), r_uhat(i));
    end
end
end
