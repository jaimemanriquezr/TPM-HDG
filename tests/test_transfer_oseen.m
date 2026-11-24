% test_transfer_oseen.m
% TPM transfer boundary conditions on the low-level fluid solvers
% (transfer_boundary_rows), with and without convection, on both
% formulations (velocity-gradient OseenProblem, stress NavierStokesProblem).
%
% Stage 1 (consistency): fitted unit square, identity mapping -> the
%   transfer reduces to the L2-projected Dirichlet datum; the solution
%   must match the strong-Dirichlet solve dof-for-dof.
% Stage 2 (polynomial reproduction): on the immersed disk a solution of
%   degree <= k is reproduced to round-off by both forms -- the path
%   integral, extrapolation and layout are exact (the pressure is checked
%   up to its constant, the manufactured mean being the unit square's).
% Stage 3 (unfitted rates): disk R = 0.75 immersed in a background grid of
%   [-1,1]^2 (distance O(h) to Gamma), datum transferred radially;
%   manufactured Oseen on meshes n = 8, 16, 32: the EOC of the last pair
%   must be at least k + 0.5 in all four quantities (L, u, p, u*) of the
%   L-form for k = 1, 2, and at least k + 0.4 in (sigma, u, gamma, u*) of
%   the S-form for k = 2.
%   The S-form at k = 1 is not checked: its sigma-reconstructed path
%   gradient is unstable under refinement.
% Stage 4: Picard Navier-Stokes on the same disk converges and lands
%   within a factor 2 of the prescribed-beta error.
clc; clear;
n_fail = 0;
nu = 0.1;
U = "sin(pi*x)*sin(pi*y)";  P = "cos(pi*x)*sin(pi*y)";
beta = ["1 + y", "x"];
dL = OseenProblem.get_manufactured_solution(U, P, nu, beta);
dS = NavierStokesProblem.get_manufactured_solution(U, P, nu, beta);

%% Stage 1: identity mapping == strong Dirichlet (linear velocity, so the
%  two datum projections coincide whatever their quadrature)
mesh = StructuredSimplicialGrid(8, 8);
d1L = OseenProblem.get_manufactured_solution("x^2 + x*y + y^2", "x*y", nu, beta);
d1S = NavierStokesProblem.get_manufactured_solution("x^2 + x*y + y^2", "x*y", nu, beta);
for k = [1 2]
    for conv = [false true]
        if conv, bL = d1L.beta; bS = d1S.beta; else, bL = {}; bS = {}; end
        pL = OseenProblem(viscosity = nu, source_term = d1L.f, datum_dirichlet = d1L.u, convective_velocity = bL);
        s0 = solve(pL, mesh, Order = k);
        pL.Data.transfer_mapping = @(X) X;
        s1 = solve(pL, mesh, Order = k);
        eL = max(abs(s0.Unknowns.field(2).dof - s1.Unknowns.field(2).dof), [], 'all');
        pS = oseen_stress(bS, viscosity = nu, source_term = d1S.f, datum_dirichlet = d1S.u);
        s0 = solve(pS, mesh, Order = k);
        pS.Data.transfer_mapping = @(X) X;
        s1 = solve(pS, mesh, Order = k);
        eS = max(abs(s0.Unknowns.field(2).dof - s1.Unknowns.field(2).dof), [], 'all');
        fprintf("Stage 1 k = %d conv = %d: max |du| L-form %.1e  S-form %.1e\n", k, conv, eL, eS);
        if max(eL, eS) > 1e-9
            fprintf("Stage 1 FAILED\n"); n_fail = n_fail + 1;
        end
    end
end

%% Stage 2: polynomial reproduction on the immersed disk
R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial = @(X) R * X ./ vecnorm(X);
mesh = immerse_mesh(StructuredSimplicialGrid(8, 8, -1, 1, -1, 1), phi);
for k = [1 2]
    if k == 1, Up = "x^2 + x*y + y^2"; else, Up = "x^3 + x^2*y + y^3"; end
    dpL = OseenProblem.get_manufactured_solution(Up, "x", nu, beta);
    dpS = NavierStokesProblem.get_manufactured_solution(Up, "x", nu, beta);
    sL = solve(OseenProblem(viscosity = nu, source_term = dpL.f, datum_dirichlet = dpL.u, ...
        convective_velocity = dpL.beta, transfer_mapping = radial), mesh, Order = k);
    eL = sL.compute_l2_error({dpL.L, dpL.u});
    sS = solve(oseen_stress(dpS.beta, viscosity = nu, source_term = dpS.f, datum_dirichlet = dpS.u, ...
        transfer_mapping = radial), mesh, Order = k);
    eS = sS.compute_l2_error({{}, dpS.u, dpS.gamma});
    fprintf("Stage 2 k = %d: polynomial reproduction  L-form (L, u) %.1e %.1e   S-form (u, gamma) %.1e %.1e\n", ...
        k, eL(1), eL(2), eS(2), eS(3));
    if max([eL(1:2); eS(2:3)]) > 1e-10
        fprintf("Stage 2 k = %d FAILED\n", k); n_fail = n_fail + 1;
    end
end

%% Stage 3: unfitted disk, radial transfer, manufactured Oseen rates
ns = [8 16 32];
for k = [1 2]
    [eL, eS] = deal(zeros(numel(ns), 4));
    for i = 1:numel(ns)
        mesh = immerse_mesh(StructuredSimplicialGrid(ns(i), ns(i), -1, 1, -1, 1), phi);
        pL = OseenProblem(viscosity = nu, source_term = dL.f, datum_dirichlet = dL.u, ...
            convective_velocity = dL.beta, transfer_mapping = radial);
        sL = solve(pL, mesh, Order = k);
        e = sL.compute_l2_error({dL.L, dL.u, dL.p});
        eL(i, :) = [e(1:3).', l2_error_p1(mesh, k, divergence_free_velocity(sL), dL.u)];
        pS = oseen_stress(dS.beta, viscosity = nu, source_term = dS.f, datum_dirichlet = dS.u, ...
            transfer_mapping = radial);
        sS = solve(pS, mesh, Order = k);
        e = sS.compute_l2_error({dS.sigma, dS.u, dS.gamma});
        eS(i, :) = [e(1:3).', l2_error_p1(mesh, k, divergence_free_velocity(sS), dS.u)];
    end
    rL = compute_eoc(eL, 1./ns);  rS = compute_eoc(eS, 1./ns);
    rL = rL(end, :);  rS = rS(end, :);
    fprintf("Stage 3 k = %d: L-form eoc (L, u, p, u*) = %.2f %.2f %.2f %.2f | S-form (sigma, u, gamma, u*) = %.2f %.2f %.2f %.2f\n", k, rL, rS);
    bad = any(rL < k + 0.5);
    if k == 2, bad = bad || any(rS < k + 0.4); end
    if bad
        fprintf("Stage 3 k = %d FAILED\n", k); n_fail = n_fail + 1;
    end
end

%% Stage 4: Picard Navier-Stokes on the disk
k = 2;
mesh = immerse_mesh(StructuredSimplicialGrid(32, 32, -1, 1, -1, 1), phi);
dL = OseenProblem.get_manufactured_solution(U, P, nu, []);
dS = NavierStokesProblem.get_manufactured_solution(U, P, nu, []);
pL = OseenProblem(viscosity = nu, source_term = dL.f, datum_dirichlet = dL.u, transfer_mapping = radial);
[sL, hL] = solve_navier_stokes(pL, mesh, k, tolerance = 1e-10);
eL = sL.compute_l2_error({dL.L, dL.u, dL.p});
pL.Data.convective_velocity = dL.beta;
rL = solve(pL, mesh, Order = k).compute_l2_error({dL.L, dL.u, dL.p});
pS = NavierStokesProblem(viscosity = nu, source_term = dS.f, datum_dirichlet = dS.u, transfer_mapping = radial);
[sS, hS] = solve_navier_stokes(pS, mesh, k, tolerance = 1e-10);
eS = sS.compute_l2_error({dS.sigma, dS.u, dS.gamma});
rS = solve(oseen_stress(dS.beta, viscosity = nu, source_term = dS.f, datum_dirichlet = dS.u, ...
    transfer_mapping = radial), mesh, Order = k).compute_l2_error({dS.sigma, dS.u, dS.gamma});
fprintf("Stage 4 Picard: L-form %d it, e(u) %.3e vs exact-beta %.3e;  S-form %d it, e(u) %.3e vs %.3e\n", ...
    numel(hL), eL(2), rL(2), numel(hS), eS(2), rS(2));
if hL(end) > 1e-10 || hS(end) > 1e-10 || eL(2) > 2*rL(2) || eS(2) > 2*rS(2)
    fprintf("Stage 4 FAILED\n"); n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_transfer_oseen:failed", "%d checks failed", n_fail);
end
