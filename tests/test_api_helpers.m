% test_api_helpers.m
% Public helpers outside the solve paths do what they declare. One check
% per helper:
%   - basis-level compute_l2_error (Dubiner, Legendre, PEERS) equals
%     FEField.compute_l2_error;
%   - the generic AbstractFEFunction.project equals Datum.project;
%   - Datum.project onto PEERS is exact on P_k tensors, rate k+1;
%   - deviator(A) = A - tr(A)/2 I;
%   - IndexedMatrix diag (block diagonal of the element pages) and
%     reshape (keeps the ids);
%   - IndexedMatrix * FEField is the local matrix-vector product;
%   - ParametrizedCurve.get_partition_area is the arc length of each patch;
%   - SimplicialGrid.get_normal_by_element is the outward normal;
%   - FluidStructureCoupling.solve_steady equals solve(..., "steady");
%   - dot(beta, n)*u_hat (either order) assembles as the Product
%     u_hat*beta.'*n on the element boundary;
%   - a test expression without data (v, c*v, bubble*v) times a measure
%     is a linear form.
clc; clear;

n_fail = 0;
mesh = StructuredSimplicialGrid(4, 4, 0, 1, 0, 1);
n_el = size(mesh.elements, 2);
f = @(x, y) sin(x).*exp(y);
k = 2;

%% basis-level compute_l2_error
B = DubinerBasis(k); L = LegendreBasis(k); P = PEERSBasis(k);
uh = Datum(f).project(mesh, B);
lh = Datum(f).project(mesh, L);
F = {@(x, y) sin(x).*y, @(x, y) x.^2; @(x, y) exp(x), @(x, y) cos(y)};
sh = FEField(P, sin(1:P.Dimension*n_el));
E_P = P.compute_l2_error(sh.dof, mesh, exact_solution = F);
e = [B.compute_l2_error(uh.dof, mesh, exact_solution = f) / uh.compute_l2_error(mesh, Datum(f)), ...
     L.compute_l2_error(lh.dof, mesh, exact_solution = f) / lh.compute_l2_error(mesh, Datum(f)), ...
     norm(E_P(:)) / sh.compute_l2_error(mesh, Datum(F))] - 1;
n_fail = n_fail + check(max(abs(e)) < 1e-12, "basis compute_l2_error = FEField.compute_l2_error");

%% generic project = Datum.project (scalar, face, vector)
VB = B.to_vector();
G = {@(x, y) sin(x).*y; @(x, y) x.^2};
e = [max(abs(project(2*Datum(f), mesh, B).dof - Datum(@(x, y) 2*f(x, y)).project(mesh, B).dof), [], "all"), ...
     max(abs(project(2*Datum(f), mesh, L).dof - Datum(@(x, y) 2*f(x, y)).project(mesh, L).dof), [], "all"), ...
     max(abs(project(2*Datum(G), mesh, VB).dof ...
        - Datum({@(x, y) 2*sin(x).*y; @(x, y) 2*x.^2}).project(mesh, VB).dof), [], "all")];
n_fail = n_fail + check(max(e) < 1e-12, "AbstractFEFunction.project = Datum.project");

%% Datum.project onto PEERS: exact on P_k tensors, rate k+1 on smooth data
[~, e_poly] = Datum({@(x, y) x.*y, @(x, y) 1 + 0*x; @(x, y) x.^2, @(x, y) y}).project(mesh, P);
Ns = [4 8 16];
e_k = zeros(3, numel(Ns));
for kk = 0:2
    for i = 1:numel(Ns)
        [~, e_k(kk + 1, i)] = Datum(F).project(StructuredSimplicialGrid(Ns(i), Ns(i)), PEERSBasis(kk));
    end
end
eoc = compute_eoc(e_k.', 1./Ns);       % one row per level, one column per k
n_fail = n_fail + check(e_poly < 1e-12 && all(abs(eoc(end, :) - (1:3)) < 0.1), ...
    "Datum.project onto PEERS (exact on P_k, rate k+1)");

%% deviator
[xq, ~] = get_quadrature(B);
a = Datum(F).evaluate(xq, mesh);
d = deviator(Datum(F)).evaluate(xq, mesh);
tr = a(:, :, 1, 1, :) + a(:, :, 2, 2, :);
n_fail = n_fail + check(max(abs(d - (a - tr .* permute(eye(2)/2, [3 4 1 2]))), [], "all") < 1e-14, ...
    "deviator(A) = A - tr(A)/2 I");

%% IndexedMatrix diag and reshape
dx = Measure(preset="dx");
u1 = TrialFunction(DubinerBasis(1)); v1 = TestFunction(DubinerBasis(1));
M = IndexedMatrix(assemble(u1*v1*dx, mesh));
m = double(M);
pages = num2cell(reshape(m, 3, 3, []), [1 2]);
R = reshape(IndexedMatrix(M.Arrays{1}, IndexedArray(2*m, [1 2])), 3, 3, []);
n_fail = n_fail + check(isequal(diag(M), blkdiag(pages{:})) && isequal(R.id, {[-1 -1]; [1 2]}) ...
    && isequal(double(R.Arrays{2}), reshape(2*m, 3, 3, [])), "IndexedMatrix diag and reshape");

%% IndexedMatrix * FEField: local matrix-vector product
dS = Measure(preset="dS");
A = IndexedMatrix(assemble(TrialFunction(LegendreBasis(1))*TestFunction(DubinerBasis(1))*dS, mesh));
a = double(A);
g = FEField(FESpace(LegendreBasis(1), mesh), cos(1:2*size(mesh.faces, 2)));
ref = zeros(3, 1, 1, 1, n_el);
for el = 1:n_el
    for s = 1:3
        ref(:, 1, 1, 1, el) = ref(:, 1, 1, 1, el) + a(:, :, 1, 1, el, 1, s)*g.dof(:, mesh.boundaries(s, el));
    end
end
n_fail = n_fail + check(max(abs(double(A*g) - ref), [], "all") < 1e-15, "IndexedMatrix * FEField");

%% get_partition_area: arc length of the patches
circle = CircularCurve([.5; .5], .25);
[fm, sm] = interface_meshes(circle, "fitted", 0.1);
pF = get_partitions(circle, fm, sm);
a = circle.get_partition_area(pF);
n_fail = n_fail + check(abs(sum(a) - pi/2) < 1e-12 && max(abs(a - pF.area)) < 1e-3, ...
    "get_partition_area = patch arc length");

%% get_normal_by_element: outward, scaled by the face length
[E, S] = ndgrid(1:n_el, 1:3);
N = mesh.get_normal_by_element(E(:), S(:));
U = reshape(permute(squeeze(mesh.exterior_unit_normals), [3 2 1]), [], 2);
n_fail = n_fail + check(max(abs(N - U .* reshape(mesh.area(mesh.boundaries).', [], 1)), [], "all") < 1e-15, ...
    "get_normal_by_element = exterior normal");

%% solve_steady = solve(..., "steady")
flat = FlatCurveHorizontal(0.5);
[fm, sm] = fitted_split(flat, [1/4 1/4], 0);
[fp, sp, sources] = FluidStructureCoupling.load_problems(1, flat, Inf, 1.0, ...
    "x", "0", ["1"; "0"], 1, 1.0, 1, 1.0, 0.25, 0.0);
[fs, ss] = FluidStructureCoupling.solve(fp, sp, fm, sm, flat, sources{:}, "steady");
[fs2, ss2] = FluidStructureCoupling.solve_steady(fp, sp, fm, sm, flat, sources{:});
n_fail = n_fail + check(isequal({fs.Unknowns.field.dof}, {fs2.Unknowns.field.dof}) ...
    && isequal({ss.Unknowns.field.dof}, {ss2.Unknowns.field.dof}), "solve_steady = solve(..., ""steady"")");

%% dot(beta, n)*u_hat = Product u_hat*beta.'*n on the element boundary
n_ = ExteriorNormal();
u_hat = TrialFunction(LegendreBasis(1).to_vector(2));
v = TestFunction(DubinerBasis(1).to_vector());
beta = Datum({@(x, y) 1 + x; @(x, y) y.^2});
ref = Product([u_hat, beta.', n_], @AbstractFEFunction.mtimes_);
ref.signature = u_hat.signature;
A_ref = double(IndexedMatrix(assemble(dot(ref, v)*dS, mesh)));
n_fail = n_fail + check(isequal(double(IndexedMatrix(assemble(dot(dot(beta, n_)*u_hat, v)*dS, mesh))), A_ref) ...
    && isequal(double(IndexedMatrix(assemble(dot(u_hat*dot(beta, n_), v)*dS, mesh))), A_ref), ...
    "dot(beta, n)*u_hat = Product u_hat*beta.'*n");

%% a test expression without data is a linear form: v*dx, c*v*dx, w*v*dx
v2 = TestFunction(B);
L = @(form) double(IndexedMatrix(assemble(form, mesh)));
[xq, wq] = get_quadrature(B);
b_ref = (double(B.evaluate(xq)).' * (wq(:) .* bubble_eval(xq(1, :), xq(2, :)).')) .* reshape(mesh.volume, 1, []);
n_fail = n_fail + check(isequal(L(v2*dx), L(Datum(@(x, y) 1 + 0*x)*v2*dx)) ...
    && isequal(L(2*v2*dx), L(Datum(@(x, y) 2 + 0*x)*v2*dx)) ...
    && max(abs(squeeze(L(BubbleFunction()*v2*dx)) - b_ref), [], "all") < 1e-15 ...
    && max(abs(L(BubbleFunction()*v2*dS)), [], "all") < 1e-15, ...
    "v*dx, c*v*dx and bubble*v*dx are linear forms");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_api_helpers:failed", "%d checks failed", n_fail);
end
