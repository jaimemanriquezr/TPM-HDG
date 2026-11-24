%% Derivatives of FEFields with dof, and their plots
% gradient(u) and divergence(v) of a field with dof evaluate to the
% field's own derivative (checked against centred differences of
% FEField.evaluate on a triangle mesh, and exactly for u in Q^2 on a
% CartesianGrid), while for trial functions (no dof) they still give one
% column per basis function. plot(gradient(u), mesh) and
% plot(divergence(v), mesh) draw one patch per element and component,
% with the field derivative at the element vertices.
clc; clear;
n_fail = 0;
rng(3);

%% triangle mesh, against centred differences
mesh = StructuredSimplicialGrid(3, 2);
n_el = size(mesh.elements, 2);
u = FEField(DubinerBasis(2), rand(6, n_el));             u.mesh = mesh;
v = FEField(DubinerBasis(2).to_vector(), rand(12, n_el)); v.mesh = mesh;
rs = [0.2 0.5; 0.3 0.1];
[gu_fd, dv_fd] = centred_differences(u, v, rs, mesh);
gu = gradient(u).evaluate(rs, mesh);                    % n_pts x 1 x 2 x 1 x n_el
dv = divergence(v).evaluate(rs, mesh);                  % n_pts x 1 x 1 x 1 x n_el
n_fail = n_fail + check("grad u, size", ~isequal(size(gu), [2 1 2 1 n_el]), 0);
n_fail = n_fail + check("grad u vs FD", max(abs(squeeze(gu) - gu_fd), [], 'all'), 1e-7);
n_fail = n_fail + check("div v, size", ~isequal(size(dv), [2 1 1 1 n_el]), 0);
n_fail = n_fail + check("div v vs FD", max(abs(squeeze(dv) - dv_fd), [], 'all'), 1e-7);
gu_default = gradient(u).evaluate(rs);                  % mesh defaults to u.mesh
n_fail = n_fail + check("default mesh", max(abs(gu_default - gu), [], 'all'), 0);

% trial functions keep one column per basis function
gt = gradient(TrialFunction(DubinerBasis(2))).evaluate(rs, mesh);
dt = divergence(TrialFunction(DubinerBasis(2).to_vector())).evaluate(rs, mesh);
n_fail = n_fail + check("trial grad columns", abs(size(gt, 2) - 6), 0);
n_fail = n_fail + check("trial div columns", abs(size(dt, 2) - 12), 0);
% ... and contracting them by hand gives the field derivative
gt_dof = sum(gt .* reshape(u.dof, 1, 6, 1, 1, []), 2);
n_fail = n_fail + check("trial . dof = field", max(abs(gt_dof - gu), [], 'all'), 1e-13);

%% CartesianGrid, exact for u in Q^2
grid = CartesianGrid([0 0.3 1], [0 0.5 0.8 1]);
f = @(x, y) x.^2 .* y + 2*x.*y.^2;
[xy, rq, w] = quadpoints_elements_2d(grid, 3, output = "array");
Pq = QkBasis(2).evaluate(rq);
uq = FEField(QkBasis(2), Pq.' * (w .* reshape(f(xy(1,:,:), xy(2,:,:)), [], 6)));  % orthonormal
uq.mesh = grid;
X = grid.get_physical_points(rs);
gx = reshape(2*X(1,:,:).*X(2,:,:) + 2*X(2,:,:).^2, size(rs, 2), 1, []);
gy = reshape(X(1,:,:).^2 + 4*X(1,:,:).*X(2,:,:), size(rs, 2), 1, []);
gq = gradient(uq).evaluate(rs, grid);                   % n_pts x 1 x 2 x 1 x n_el
n_fail = n_fail + check("Q^2 grad exact", ...
    max(abs(squeeze(gq) - [gx, gy]), [], 'all'), 1e-12);

%% plots: one patch per element and component, vertex values = derivative
cases = {gradient(u), mesh, 2, 3; divergence(v), mesh, 1, 3; gradient(uq), grid, 2, 4};
for c = 1:size(cases, 1)
    [fun, m, n_comp, n_vert] = cases{c, :};
    fig = figure(Visible = "off");
    plot(fun, m);
    p = findobj(fig, "Type", "patch");
    vals = fun.evaluate(reference_element(2, m.ElementType), m);
    stored = sort(cell2mat(arrayfun(@(h) h.Vertices(:, 3), p, UniformOutput = false)));
    n_fail = n_fail + check(sprintf("plot %d patches", c), ...
        abs(numel(p) - n_comp*size(m.elements, 2)) + any(arrayfun(@(h) size(h.Vertices, 1), p) ~= n_vert), 0);
    n_fail = n_fail + check(sprintf("plot %d values", c), max(abs(stored - sort(vals(:)))), 1e-13);
    close(fig);
end
fprintf("%d check(s) failed\n", n_fail);

function [gu, dv] = centred_differences(u, v, rs, mesh)
% physical gradient of u and divergence of v at the reference points rs of
% every element, by centred differences of FEField.evaluate
h = 1e-6;
G = mesh.get_affine_map(inverse = true);
[n_pts, n_el] = deal(size(rs, 2), size(mesh.elements, 2));
[gu, dv] = deal(zeros(n_pts, 2, n_el), zeros(n_pts, n_el));
for K = 1:n_el
    for c = 1:2
        e = G(:, c, K);                                  % d(rs)/dx_c
        dU = (u.evaluate(rs + h*e) - u.evaluate(rs - h*e)) / (2*h);
        dV = (v.evaluate(rs + h*e) - v.evaluate(rs - h*e)) / (2*h);
        gu(:, c, K) = dU(:, 1, 1, 1, K);
        dv(:, K) = dv(:, K) + dV(:, 1, c, 1, K);
    end
end
end

function failed = check(name, err, tol)
    failed = ~(err <= tol);
    fprintf("  %-22s err = %.2e  %s\n", name, err, string(~failed));
end
