% test_darcy_convergence_check.m
% Darcy HDG convergence check: L2 error of the scalar unknown p on
% successively refined meshes; expected EOC k+1.
clc; clear;
n_fail = 0;

k = 2;
tau = 1;
kappa = eye(2);

filename = cache_dir() + "/darcy_solutions";
data = DarcyProblem.get_manufactured_solution("exp(sin(pi*x*y))", kappa, filename + ".txt");

n_x = @(x, y) sign(x - 1/2) .* (abs(x - 1/2) == 1/2);
n_y = @(x, y) sign(y - 1/2) .* (abs(y - 1/2) == 1/2);
flux_bc = @(u) @(x, y) u{1}(x, y).*n_x(x, y) + u{2}(x, y).*n_y(x, y);

problem = DarcyProblem( ...
    diffusion       = kappa, ...
    source          = data.f, ...
    datum_dirichlet = data.p, ...
    datum_neumann   = flux_bc(data.u));

n_ref = 3;
[err_p, h_list] = deal(zeros(n_ref, 1));
for i = 1:n_ref
    n = 4 * 2^(i-1);
    h_list(i) = 1/n;
    mesh = StructuredSimplicialGrid(n, n).add_stabilization(tau);
    solution = solve(problem, mesh, Order= k);

    err_p(i) = compute_l2_error(solution.Unknowns.field(2), mesh, Datum(data.p));
    fprintf("n = %2d, h = %.4f, L2 err p = %.4e\n", n, h_list(i), err_p(i));
end
eoc = compute_eoc(err_p, h_list);
fprintf("EOC(p): %s   (expected ~ %d)\n", num2str(eoc.', "%.2f "), k+1);

% rate between the coarsest and finest levels (the mean EOC)
rate = log(err_p(1)/err_p(end)) / log(h_list(1)/h_list(end));
fprintf("rate(p) over all levels: %.2f (pass >= %.2f)\n", rate, k + 1 - 0.25);
if rate < k + 1 - 0.25
    fprintf("rate(p) FAILED\n"); n_fail = n_fail + 1;
end

fprintf("\n%d check(s) failed.\n", n_fail);
