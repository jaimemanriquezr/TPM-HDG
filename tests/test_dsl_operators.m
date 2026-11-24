% test_dsl_operators.m
% Operator semantics of the form DSL: products are evaluated in the
% written order, while the trial/test roles come from the signatures.
%
% Check 1: test-left scalar products equal their dot spelling.
% Check 2: roles do not depend on placement (dot, inner, sum(.*), u'*v).
% Check 3: non-1x1 values in the written order (T*u, u'*T', v*u.').
% Check 4: inner-dimension mismatch errors with BilinearForm:mtimes.
% Check 5: sums of trial functions split per summand in either placement.
% Check 6: sums of test functions split per summand (bilinear and linear).
% Check 7: data contracts a non-1x1 form after the fact (re-association).
% Check 8: function-handle products resolve independently of the operand.
% Check 9: Darcy solve with the sum spellings (p - p_hat)*w, w*(p - p_hat).
% Check 10: native vs dot spelling of the FSI fluid block, term by term.
clc; clear;

n_fail = 0;
check = @(name, ok) fprintf("%-62s %s\n", name, string(ok));

k = 1; mesh = StructuredSimplicialGrid(2, 2).add_stabilization(1);
Pk = DubinerBasis(k);
dx = Measure(preset = "dx"); dS = Measure(preset = "dS");
n_ = ExteriorNormal(); tau_ = Stabilization();
S = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
                        FESpace(LegendreBasis(k).to_vector(2))]);
[sg, u, xi, u_hat] = get_trial_functions(S.Space);
[T, v, z, v_hat]   = get_test_functions(S.Space);
A   = @(F, m) double(assemble(F*m, mesh));
rel = @(X, Y) max(abs(X - Y), [], 'all') / max(max(abs(Y), [], 'all'), eps);
% blockwise comparison of two assembled IndexedMatrix on the same id set
blk = @(Y, id) double(Y.Arrays{ismember(cell2mat(Y.id), id, 'rows')});
same_ids = @(Y, Z) isequal(sortrows(cell2mat(Y.id)), sortrows(cell2mat(Z.id)));
cmp = @(Y, Z) max(arrayfun(@(i) rel(blk(Y, Y.id{i}), blk(Z, Y.id{i})), 1:numel(Y.id)));
fv = Datum({@(x, y) x.*y; @(x, y) x + y});

%% 1. test-left scalar products
r = [rel(A(v'*u, dx), A(dot(u, v), dx)), ...
     rel(A(divergence(T)'*u, dx), A(dot(u, divergence(T)), dx)), ...
     rel(A(v_hat'*(sg*n_), dS), A(dot(sg*n_, v_hat), dS)), ...
     rel(A(trace(T)*xi, dx), A(xi*trace(T), dx)), ...
     rel(A(v'*(tau_*u), dS), A(dot(tau_*u, v), dS))];
ok = all(r <= 1e-14); check("1. test-left products == dot spelling", ok);
n_fail = n_fail + ~ok;

%% 2. roles from signatures
r = [rel(A(dot(v, u), dx), A(dot(u, v), dx)), rel(A(dot(T, sg), dx), A(dot(sg, T), dx)), ...
     rel(A(sum(sg.*T), dx), A(sum(T.*sg), dx)), rel(A(sum(T.*sg), dx), A(dot(sg, T), dx)), ...
     rel(A(u'*v, dx), A(v'*u, dx)), ...
     rel(A(inner(u, v), dx), A(dot(u, v), dx)), rel(A(inner(v, u), dx), A(dot(u, v), dx))];
ok = all(r == 0); check("2. placement-independent roles (dot, inner, sum, u'*v)", ok);
n_fail = n_fail + ~ok;

%% 3. non-1x1 values in the written order
X = A(T*u, dx);
ok = isequal(size(X, 3:4), [2 1]);
for c = 1:2
    ok = ok && rel(X(:, :, c, 1, :), A((double([c == 1, c == 2])*T)*u, dx)) <= 1e-13;
end
Y = A(u'*T', dx);
ok = ok && isequal(size(Y, 3:4), [1 2]) && rel(permute(Y, [1 2 4 3 5]), X) <= 1e-13;
O = A(v*u.', dx);
ok = ok && isequal(size(O, 3:4), [2 2]) && rel(O, permute(A(u*v.', dx), [1 2 4 3 5])) == 0;
check("3. T*u, u'*T', v*u.' in the written order", ok);
n_fail = n_fail + ~ok;

%% 4. inner-dimension mismatch
try
    assemble(v*u*dx, mesh);
    ok = false;
catch ME
    ok = ME.identifier == "BilinearForm:mtimes";
end
check("4. v*u errors with BilinearForm:mtimes", ok);
n_fail = n_fail + ~ok;

%% 5. trial sums, either placement
Ya = assemble(v'*(u - u_hat)*dS, mesh);
Yb = assemble(v'*u*dS - v'*u_hat*dS, mesh);
Yc = assemble((u - u_hat)'*v*dS, mesh);
Yd = assemble((u - u_hat)*v'*dS, mesh);
ok = same_ids(Ya, Yb) && same_ids(Yc, Yb) && cmp(Ya, Yb) <= 1e-14 && cmp(Yc, Yb) <= 1e-14 ...
    && isequal(sortrows(cell2mat(Ya.id)), [2 2; 4 2]) ...
    && all(cellfun(@(Z) isequal(size(Z, 3:4), [2 2]), Yd.Arrays));
check("5. v'*(u - u_hat), (u - u_hat)'*v split per summand", ok);
n_fail = n_fail + ~ok;

%% 6. test sums
Ta = assemble((v - v_hat)'*u*dS, mesh);
Tb = assemble(v'*u*dS - v_hat'*u*dS, mesh);
Tc = assemble((v - v_hat)'*(u - u_hat)*dS, mesh);
Td = assemble(v'*u*dS - v'*u_hat*dS - v_hat'*u*dS + v_hat'*u_hat*dS, mesh);
La = assemble((v - v_hat)'*fv*dS, mesh);
Lb = assemble(v'*fv*dS - v_hat'*fv*dS, mesh);
Lc = assemble(fv'*(v - v_hat)*dS, mesh);
ok = same_ids(Ta, Tb) && cmp(Ta, Tb) <= 1e-14 ...
    && same_ids(Tc, Td) && size(cell2mat(Tc.id), 1) == 4 && cmp(Tc, Td) <= 1e-14 ...
    && same_ids(La, Lb) && cmp(La, Lb) <= 1e-14 && cmp(Lc, Lb) <= 1e-14;
check("6. (v - v_hat)'*u, (v - v_hat)'*(u - u_hat), (v - v_hat)'*f split", ok);
n_fail = n_fail + ~ok;

%% 7. contraction by data == re-associated spelling
ref = A((T'*n_)'*u, dS);
e1 = A(([1 0]*T)*u, dx);
r = [rel(A(n_'*(T*u), dS), ref), rel(A((T*u)'*n_, dS), ref), ...
     rel(A(n_.'*(T*u), dS), ref), ...
     rel(A([1 0]*(T*u), dx), e1), rel(A((u'*T')*[1; 0], dx), e1)];
ok = all(r <= 1e-14) && norm(ref(:)) > 0;
check("7. n_'*(T*u), (T*u)'*n_, [1 0]*(T*u), (u'*T')*e1 contract", ok);
n_fail = n_fail + ~ok;

%% 8. function-handle products: same route for trial and test operands
f = @(x, y) 1 + x.^2 + 2*y;
ref_u = A(dot(u*Datum(f), v), dx);
ref_z = A(z*Datum(f), dx);
r = [rel(A(dot(u*f, v), dx), ref_u), rel(A(dot(f*u, v), dx), ref_u), ...
     rel(A(v'*(u*f), dx), ref_u), ...
     rel(A(z*f, dx), ref_z), rel(A(f*z, dx), ref_z), ...
     rel(A((xi*f)*z, dx), A(z*(f*xi), dx))];
ok = all(r <= 1e-14) && norm(ref_u(:)) > 0;
check("8. u*@f, @f*u, w*@f, @f*w resolve alike", ok);
n_fail = n_fail + ~ok;

%% 9. Darcy with sums of trial functions (as tests/test_dsl_darcy.m)
kd = 2; Nd = 8;
pfun = @(x, y) sin(pi*x).*sin(pi*y);
p_ = Datum(pfun); f_ = Datum(@(x, y) 2*pi^2*sin(pi*x).*sin(pi*y));
Pd = DubinerBasis(kd);
mesh_d = StructuredSimplicialGrid(Nd, Nd).add_stabilization(1);
problem = FEProblem(Spaces = [FESpace([Pd.to_vector(), Pd], Condense = true), ...
                              FESpace(LegendreBasis(kd))]);
problem.BoundaryConditions = WeakEssentialBoundaryConditions(p_, "Dirichlet");
[q, p, p_hat] = get_trial_functions(problem.Space);
[r_, w, w_hat] = get_test_functions(problem.Space);
n = ExteriorNormal();
flux = -(dot(q, n)*w_hat*dS + p*w_hat*dS - p_hat*w_hat*dS);
local = @(jump) dot(q, r_)*dx - p*divergence(r_)*dx + p_hat*dot(r_, n)*dS ...
    + divergence(q)*w*dx + jump + f_*w*dx;
problem.Operator = AffineOperator(local(p*w*dS - p_hat*w*dS)) + flux;
sol_ref = solve(problem, mesh_d, kd);
problem.Operator = AffineOperator(local((p - p_hat)*w*dS)) + flux;
sol_a = solve(problem, mesh_d, kd);
problem.Operator = AffineOperator(local(w*(p - p_hat)*dS)) + flux;
sol_b = solve(problem, mesh_d, kd);
d = @(s) max(abs(s.Unknowns.field(2).dof - sol_ref.Unknowns.field(2).dof), [], 'all');
ok = d(sol_a) < 1e-10 && d(sol_b) < 1e-10;
check(sprintf("9. Darcy (p - p_hat)*w, w*(p - p_hat): %.1e, %.1e", d(sol_a), d(sol_b)), ok);
n_fail = n_fail + ~ok;

%% 10. FSI fluid block: native vs dot spelling, term by term
J = [0, 1; -1, 0]; mu = 0.7; c1 = 1/(2*mu); c2 = 1/(4*mu);
F = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
                        FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
[sg, u, gm, xi, u_hat, rho] = get_trial_functions(F.Space);
[T, v, e, z, v_hat, s_rho]  = get_test_functions(F.Space);
native = {c1*sum(T.*sg)*dx, c2*trace(T)*trace(sg)*dx, divergence(T)'*u*dx, ...
    sum(T.*(gm*J))*dx, trace(T)*xi*dx, sum(T.*(u_hat*n_.'))*dS, ...
    v'*divergence(sg)*dx, mu*v'*(tau_*u)*dS, mu*v'*(tau_*u_hat)*dS, v'*fv*dx, ...
    sum((e*J).*sg)*dx, z*trace(sg)*dx, z*rho*dx, ...
    v_hat'*(sg*n_)*dS, mu*v_hat'*(tau_*u)*dS, mu*v_hat'*(tau_*u_hat)*dS, s_rho*(u_hat'*n_)*dS};
dotted = {dot(sg, c1*T)*dx, trace(sg)*(c2*trace(T))*dx, dot(u, divergence(T))*dx, ...
    dot(gm*J, T)*dx, xi*trace(T)*dx, dot(u_hat*n_.', T)*dS, ...
    dot(divergence(sg), v)*dx, dot((mu*tau_)*u, v)*dS, dot((mu*tau_)*u_hat, v)*dS, dot(fv, v)*dx, ...
    dot(sg, e*J)*dx, trace(sg)*z*dx, rho*z*dx, ...
    dot(sg*n_, v_hat)*dS, dot((mu*tau_)*u, v_hat)*dS, dot((mu*tau_)*u_hat, v_hat)*dS, dot(u_hat, n_)*s_rho*dS};
r = cellfun(@(a, b) rel(double(assemble(a, mesh)), double(assemble(b, mesh))), native, dotted);
ok = all(r <= 1e-14);
check(sprintf("10. FSI fluid native vs dot, %d terms (max rel %.1e)", numel(r), max(r)), ok);
n_fail = n_fail + ~ok;

fprintf("\n%d check(s) failed.\n", n_fail);
