% test_weak_bc_elasticity.m
% Baseline for the declarative solve path with the PEERS elasticity system:
% fitted unit square, weakly imposed Dirichlet data, manufactured solution.
clc; clear;
n_fail = 0;

k = 2;
E = 1.0;
nu = 0.3;
lambda = E*nu/((1 + nu)*(1 - 2*nu));
mu = E/2/(1 + nu);

syms x y
u_sym = [sin(pi*x) .* sin(pi*y); cos(pi*x) .* cos(pi*y)];
gu = [diff(u_sym(1), x), diff(u_sym(1), y);
      diff(u_sym(2), x), diff(u_sym(2), y)];
eps_u = (gu + gu.')/2;
sig = lambda*trace(eps_u)*eye(2) + 2*mu*eps_u;
f_sym = -[diff(sig(1,1), x) + diff(sig(1,2), y);
          diff(sig(2,1), x) + diff(sig(2,2), y)];

u_ = Datum(u_sym);
f1 = Datum(f_sym(1));
f2 = Datum(f_sym(2));

S_k = PEERSBasis(k);
Pk = DubinerBasis(k);
VPk = Pk.to_vector();
VMk = LegendreBasis(k).to_vector(2);

dx = Measure(preset="dx");

n_ref = 3;
[err_u, h_list] = deal(zeros(n_ref, 1));
for i = 1:n_ref
    n = 4 * 2^(i-1);
    h_list(i) = 1/n;
    mesh = StructuredSimplicialGrid(n, n);

    problem = FEProblem(Spaces= [FESpace([S_k, VPk, Pk], Condense=true), FESpace(VMk)]);
    [lhs_op, flux_op] = LinearElasticityMatrices(k, mesh, E, nu);
problem.Operator = AffineOperator(lhs_op, flux_op);

    % source moments against the vector Dubiner test space ([x; y] blocks)
    w = TestFunction(Pk);
    F1 = double(f1 * w * Measure(mesh, 'dx'));
    F2 = double(f2 * w * Measure(mesh, 'dx'));
    F = IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));
    problem.Operator = problem.Operator + F;

    problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");
    solution = solve(problem, mesh, k);

    u_h = solution.Unknowns.field(2);
    err_u(i) = compute_l2_error(u_h, mesh, u_);
    fprintf("n = %2d, L2 err u = %.4e\n", n, err_u(i));
end
eoc = compute_eoc(err_u, h_list);
fprintf("EOC(u): %s (expected ~ %d)\n", num2str(eoc.', "%.2f "), k+1);

% rate between the coarsest and finest levels (the mean EOC)
rate = log(err_u(1)/err_u(end)) / log(h_list(1)/h_list(end));
fprintf("rate(u) over all levels: %.2f (pass >= %.2f)\n", rate, k + 1 - 0.25);
if rate < k + 1 - 0.25
    fprintf("rate(u) FAILED\n"); n_fail = n_fail + 1;
end

fprintf("\n%d check(s) failed.\n", n_fail);
