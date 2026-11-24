% test_oseen.m
% Oseen / Navier–Stokes on both fluid formulations (Solano & Vargas 2022,
% eq. 3.4): (1) with beta = {} the Oseen classes reproduce their Stokes
% parents dof-for-dof; (2) a manufactured Oseen solution with prescribed
% beta converges at k+1 in L, u, p (velocity-gradient form) and in sigma,
% u (stress form); (3) the divergence-free postprocessing u* converges at
% k+2 with div u* = 0; (4) Picard on the manufactured Navier–Stokes problem
% (beta = u) converges and lands on the prescribed-beta error.
clc; clear;
n_fail = 0;

U = "sin(pi*x)^2*sin(pi*y)^2";  P = "cos(pi*x)*sin(pi*y)";  nu = 0.1;
f = {@(x, y) sin(pi*x).*cos(pi*y); @(x, y) exp(x).*y - exp(1)/2 + 1/2};

%% (1) beta = {} is Stokes
mesh = StructuredSimplicialGrid(6, 5);
for k = [1 2]
    [~, H1, R1] = get_matrices(OseenProblem(k, viscosity = nu, source_term = f), mesh);
    [~, H0, R0] = get_matrices(StokesProblem(k, viscosity = nu, source_term = f), mesh);
    [~, H3, R3] = get_matrices(NavierStokesProblem(k, viscosity = nu, source_term = f), mesh);
    [~, H2, R2] = get_matrices(IncompressibleElasticityProblem(k, viscosity = nu, source_term = f), mesh);
    if ~(isequal(H1, H0) && isequal(R1, R0) && isequal(H3, H2) && isequal(R3, R2))
        fprintf("k = %d: beta = {} does not reproduce Stokes FAILED\n", k); n_fail = n_fail + 1;
    end
end

%% (2) manufactured Oseen, prescribed beta, k+1
beta = ["1 + y", "x"];
dL = OseenProblem.get_manufactured_solution(U, P, nu, beta);
dS = NavierStokesProblem.get_manufactured_solution(U, P, nu, beta);
Ns = [8 16];
for k = [1 2]
    [eL, eS] = deal(zeros(2, 3));
    for i = 1:2
        mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
        pL = OseenProblem(viscosity = nu, source_term = dL.f, datum_dirichlet = dL.u, convective_velocity = dL.beta);
        e = solve(pL, mesh, Order = k).compute_l2_error({dL.L, dL.u, dL.p});  eL(i, :) = e(1:3);
        pS = oseen_stress(dS.beta, viscosity = nu, source_term = dS.f, datum_dirichlet = dS.u);
        e = solve(pS, mesh, Order = k).compute_l2_error({dS.sigma, dS.u, dS.gamma});  eS(i, :) = e(1:3);
    end
    rL = compute_eoc(eL, 1./Ns);  rS = compute_eoc(eS, 1./Ns);
    fprintf("k = %d: Oseen eoc  L-form (L, u, p) = %.2f %.2f %.2f   S-form (sigma, u, gamma) = %.2f %.2f %.2f\n", k, rL, rS);
    if any([rL, rS(1:2)] < k + 0.8)
        fprintf("k = %d: Oseen convergence FAILED\n", k); n_fail = n_fail + 1;
    end
end

%% (3) divergence-free postprocessing: u* at k+2, div u* = 0, both forms
k = 2;
d0L = StokesProblem.get_manufactured_solution(U, P, nu);
d0S = IncompressibleElasticityProblem.get_manufactured_solution(U, P, nu);
[eL, eS, dv] = deal(zeros(2, 1));
for i = 1:2
    mesh = StructuredSimplicialGrid(Ns(i), Ns(i));
    sL = solve(StokesProblem(viscosity = nu, source_term = d0L.f, datum_dirichlet = d0L.u), mesh, Order = k);
    [usL, dvL] = divergence_free_velocity(sL);  eL(i) = l2_error_p1(mesh, k, usL, d0L.u);
    sS = solve(IncompressibleElasticityProblem(viscosity = nu, source_term = d0S.f, datum_dirichlet = d0S.u), mesh, Order = k);
    [usS, dvS] = divergence_free_velocity(sS);  eS(i) = l2_error_p1(mesh, k, usS, d0S.u);
    dv(i) = max(dvL, dvS);
end
rL = compute_eoc(eL, 1./Ns);  rS = compute_eoc(eS, 1./Ns);
fprintf("k = %d: u* eoc L-form %.2f  S-form %.2f  (max div %.1e)\n", k, rL, rS, max(dv));
if rL < k + 1.8 || rS < k + 1.8 || max(dv) > 1e-10
    fprintf("k = %d: postprocessing FAILED\n", k); n_fail = n_fail + 1;
end

%% (4) Picard Navier–Stokes (beta = u*) lands on the prescribed-beta error
mesh = StructuredSimplicialGrid(16, 16);
dL = OseenProblem.get_manufactured_solution(U, P, nu, []);
dS = NavierStokesProblem.get_manufactured_solution(U, P, nu, []);
pL = OseenProblem(viscosity = nu, source_term = dL.f, datum_dirichlet = dL.u);
[sL, hL] = solve_navier_stokes(pL, mesh, k, tolerance = 1e-10);
eL = sL.compute_l2_error({dL.L, dL.u, dL.p});
pL.Data.convective_velocity = dL.beta;
rL = solve(pL, mesh, Order = k).compute_l2_error({dL.L, dL.u, dL.p});
pS = NavierStokesProblem(viscosity = nu, source_term = dS.f, datum_dirichlet = dS.u);
[sS, hS] = solve_navier_stokes(pS, mesh, k, tolerance = 1e-10);
eS = sS.compute_l2_error({dS.sigma, dS.u, dS.gamma});
rS = solve(oseen_stress(dS.beta, viscosity = nu, source_term = dS.f, datum_dirichlet = dS.u), ...
    mesh, Order = k).compute_l2_error({dS.sigma, dS.u, dS.gamma});
fprintf("Picard: L-form %d it, e(u) %.3e vs exact-beta %.3e;  S-form %d it, e(u) %.3e vs %.3e\n", ...
    numel(hL), eL(2), rL(2), numel(hS), eS(2), rS(2));
if hL(end) > 1e-10 || hS(end) > 1e-10 || abs(eL(2) - rL(2)) > 0.05*rL(2) || abs(eS(2) - rS(2)) > 0.05*rS(2)
    fprintf("Picard FAILED\n"); n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_oseen:failed", "%d checks failed", n_fail);
end
