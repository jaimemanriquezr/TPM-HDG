% test_dsl_darcy.m
% The full Darcy HDG system (local problem and flux equation) written
% symbolically in the form DSL must match the DarcyMatrices solution to
% machine precision.
%
% Also asserts the operand-order and spelling equivalences:
%   dot(q, n)*v_hat*dS == v_hat*dot(q, n)*dS     (Product absorbs Measure)
%   q'*n               == dot(q, n)              (ctranspose chain)
% and exercises unary minus on form arrays: flux = -(s1 + s2 - s3).
clc; clear;

k = 2; N = 16;
kappa = eye(2);
pfun = @(x, y) sin(pi*x).*sin(pi*y);
ffun = @(x, y) 2*pi^2*sin(pi*x).*sin(pi*y);
p_ = Datum(pfun); f_ = Datum(ffun);

Pk = DubinerBasis(k); VPk = Pk.to_vector(); Mk = LegendreBasis(k);
mesh = StructuredSimplicialGrid(N, N).add_stabilization(1);

n_fail = 0;

%% reference: DarcyMatrices path
problem = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)]);
problem.BoundaryConditions = WeakEssentialBoundaryConditions(p_, "Dirichlet");
[q, u, u_hat] = get_trial_functions(problem.Space);
[r, w, v_hat] = get_test_functions(problem.Space);
dx = Measure(preset="dx");
dS = Measure(preset="dS");
n = ExteriorNormal();

F = assemble(f_ * w * dx, mesh);
[lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
problem.Operator = AffineOperator(lhs_op, flux_op);
problem.Operator = problem.Operator + F;
sol_ref = solve(problem, mesh, k);

%% spelling equivalences
A1 = double(assemble(v_hat*dot(q, n)*dS, mesh));
A2 = double(assemble(dot(q, n)*v_hat*dS, mesh));
err = max(abs(A1(:) - A2(:)));
fprintf("operand order: |v_hat*dot(q,n) - dot(q,n)*v_hat|          = %.3e\n", err);
if ~(err == 0); n_fail = n_fail + 1; end

A3 = double(assemble(v_hat*(q'*n)*dS, mesh));
A4 = double(assemble((q'*n)*v_hat*dS, mesh));
err = max(abs(A1(:) - A3(:)));
fprintf("spelling:      |v_hat*dot(q,n) - v_hat*(q'*n)|            = %.3e\n", err);
if ~(err < 1e-14); n_fail = n_fail + 1; end
err = max(abs(A3(:) - A4(:)));
fprintf("operand order: |v_hat*(q'*n) - (q'*n)*v_hat|              = %.3e\n", err);
if ~(err == 0); n_fail = n_fail + 1; end

B1 = double(assemble(u_hat*dot(r, n)*dS, mesh));
B2 = double(assemble(u_hat*(r'*n)*dS, mesh));
err = max(abs(B1(:) - B2(:)));
fprintf("test side:     |u_hat*dot(r,n) - u_hat*(r'*n)|            = %.3e\n", err);
if ~(err < 1e-14); n_fail = n_fail + 1; end

%% full symbolic system (natural spelling, mixed sums, unary minus)
problem.Operator = AffineOperator(...
      dot(q, r)*dx - u*divergence(r)*dx + u_hat*dot(r, n)*dS ...
    + divergence(q)*w*dx + u*w*dS - u_hat*w*dS ...
    + f_*w*dx);
problem.Operator = problem.Operator + (-(dot(q, n)*v_hat*dS + u*v_hat*dS - u_hat*v_hat*dS));

sol_dsl = solve(problem, mesh, k);
diff_u = max(abs(sol_ref.Unknowns.field(2).dof - sol_dsl.Unknowns.field(2).dof), [], 'all');
diff_q = max(abs(sol_ref.Unknowns.field(1).dof - sol_dsl.Unknowns.field(1).dof), [], 'all');
fprintf("DSL vs DarcyMatrices: max |u dof diff| = %.3e, max |q dof diff| = %.3e\n", ...
    diff_u, diff_q);
if ~(diff_u < 1e-10 && diff_q < 1e-10); n_fail = n_fail + 1; end

err_u = compute_l2_error(sol_dsl.Unknowns.field(2), mesh, p_);
fprintf("DSL solve: L2 err p = %.6e\n", err_u);
if ~(err_u < 1e-3); n_fail = n_fail + 1; end   % k = 2, h = 1/16: ~8e-5

%% same system spelled with q'*n / r'*n
problem.Operator = AffineOperator(...
      dot(q, r)*dx - u*divergence(r)*dx + u_hat*(r'*n)*dS ...
    + divergence(q)*w*dx + u*w*dS - u_hat*w*dS ...
    + f_*w*dx);
problem.Operator = problem.Operator + (-((q'*n)*v_hat*dS + u*v_hat*dS - u_hat*v_hat*dS));

sol_ct = solve(problem, mesh, k);
diff_u = max(abs(sol_ref.Unknowns.field(2).dof - sol_ct.Unknowns.field(2).dof), [], 'all');
fprintf("q'*n spelling vs DarcyMatrices: max |u dof diff| = %.3e\n", diff_u);
if ~(diff_u < 1e-10); n_fail = n_fail + 1; end

fprintf("\n%d check(s) failed.\n", n_fail);
