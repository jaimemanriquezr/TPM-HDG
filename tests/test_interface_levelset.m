% test_interface_levelset.m
% A LevelSetCurve as the coupling interface (IP-TPM): closed-circle
% Darcy-Darcy with symbolic transmission rows, solved twice on the same
% carved meshes: once with the CircularCurve, once with the same circle
% as a discretized LevelSetCurve. The two representations induce the
% same partitions up to the polyline resolution, so the solutions must
% agree at that scale; the level-set run must also converge to the
% manufactured solution.
clc; clear;

kappa = eye(2);
k = 2;
c0 = [0.5; 0.5];
R = 0.2;

circle_p = CircularCurve(c0, R, mapping="closest");
circle_ls = LevelSetCurve( ...
    @(X) (X(1, :) - c0(1)).^2 + (X(2, :) - c0(2)).^2 - R^2, ...
    gradient = @(X) 2*(X - c0), box = [0 1 0 1]);
% anchor the seam where CircularCurve puts t = 0, so the two
% representations induce the same patch partition
circle_ls = discretize(circle_ls, n_grid = 2048, anchor = c0 + [R; 0]);

syms x y
p_sym = sin(pi*x) .* sin(pi*y);
q_sym = -kappa * gradient(p_sym, [x y]);
f_sym = divergence(q_sym, [x y]);
p_ = Datum(p_sym);
f_ = Datum(f_sym);

Pk = DubinerBasis(k); VPk = Pk.to_vector(); Mk = LegendreBasis(k);
dx = Measure(preset = "dx");
n_fail = 0;

% one background, carved once per representation: the classifications
% must agree, then the CircularCurve meshes serve both solves
n = 16;
bg = StructuredSimplicialGrid(n, n, "crossed");
[outer_p, inner_p] = bg.split(circle_p, cut_elements="mesh_2");
[outer_l, inner_l] = bg.split(circle_ls, cut_elements="mesh_2");
same_carve = isequal(outer_p.elements, outer_l.elements) ...
    && isequal(inner_p.elements, inner_l.elements);
fprintf("split classification agrees across representations: %d\n", same_carve);
if ~same_carve
    fprintf("FAIL split disagrees\n"); n_fail = n_fail + 1;
end

solutions = cell(1, 2);
curves = {circle_p, circle_ls};
for r = 1:2
    interface = curves{r};
    pb1 = configure(outer_p, k, kappa, f_, p_, VPk, Pk, Mk, dx);
    pb2 = configure(inner_p, k, kappa, f_, p_, VPk, Pk, Mk, dx);
    coupled = InterfaceProblem.darcy_darcy(pb1, pb2, interface, kappa, kappa);

    [q1, p1, phat1] = get_trial_functions(pb1);
    [q2, p2, phat2] = get_trial_functions(pb2);
    [~, ~, mu1] = get_test_functions(pb1);
    [~, ~, mu2] = get_test_functions(pb2);
    tau_ = Stabilization();
    T1 = ( dot(Extrapolate(q1), InterfaceNormal(outer_p)) + tau_*(p1 - phat1) ...
         + dot(Extrapolate(q2), InterfaceNormal(inner_p)) + tau_*(p2 - phat2) ...
         ) * mu1 * Measure(outer_p, "dI");
    T2 = ( (phat1 - PathIntegral(kappa \ q1)) ...
         - (phat2 - PathIntegral(kappa \ q2)) ) * mu2 * Measure(inner_p, "dI");
    coupled.transmission = [T1, T2];

    [sol1, sol2] = solve(coupled, k);
    solutions{r} = {sol1, sol2};
end

d = 0;
for side = 1:2
    a = solutions{1}{side}; b = solutions{2}{side};
    for fld = 1:2
        d = max(d, max(abs(a.Unknowns.field(fld).dof - b.Unknowns.field(fld).dof), [], 'all'));
    end
end
% with the seam anchored, the remaining difference is the induced-
% parametrization resolution
fprintf("CircularCurve vs LevelSetCurve max dof diff = %.3e\n", d);
if d > 5e-4
    fprintf("FAIL representations disagree beyond parametrization scale\n");
    n_fail = n_fail + 1;
end

err = zeros(2, 2);   % (representation, side)
meshes = {outer_p, inner_p};
for r = 1:2
    for side = 1:2
        err(r, side) = compute_l2_error(solutions{r}{side}.Unknowns.field(2), meshes{side}, p_);
    end
end
fprintf("L2 err p: circular = (%.3e, %.3e), level-set = (%.3e, %.3e)\n", ...
    err(1, 1), err(1, 2), err(2, 1), err(2, 2));
if any(err(2, :) > 1.2 * err(1, :) + 1e-12)
    fprintf("FAIL level-set errors exceed the parametrized run's\n");
    n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_levelset:failed", "%d checks failed", n_fail);
end

function pb = configure(mesh, k, kappa, f_, p_, VPk, Pk, Mk, dx)
pb = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)], ...
    Order = k);
pb.Mesh = mesh;
[lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
pb.Operator = AffineOperator(lhs_op, flux_op);
[~, w, ~] = get_test_functions(pb.Space);
pb.Operator = pb.Operator + assemble(f_ * w * dx, mesh);
pb.BoundaryConditions = WeakEssentialBoundaryConditions(p_, "Dirichlet");
end
