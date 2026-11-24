%-------------------------------------------------------------------------
%Description: The self-convergence error routine.
%
%             The decisive check is exactness: a polynomial of degree <= k
%             lies in both the coarse and the fine space, so the measured
%             distance between its two representations must be zero. That
%             single number exercises the parent map, the composed child
%             map and the quadrature at once, and a misalignment in any of
%             them shows up immediately -- a wrong parent gives O(1), a
%             transposed child map gives O(h^k).
%
%             The second check is that the routine reproduces a known
%             answer: against a manufactured field it must agree with
%             FEField.compute_l2_error on the same mesh.
%-------------------------------------------------------------------------
n_failed = 0;
k = 2;

coarse = StructuredSimplicialGrid(2, 2);
[fine1, info1] = uniform_refine(coarse);
[fine2, info2] = uniform_refine(coarse, 2);

%% Exactness: the same polynomial on both meshes --------------------------
% Degree k, so it is represented exactly in both spaces.
p = @(x, y) 1 - 2*x + 0.5*y + 3*x.^2 - x.*y + 0.25*y.^2;
basis = DubinerBasis(k);

u_coarse = FEField(basis, L2proj2D(coarse, p, k));
for pair = {{fine1, info1, "1 level"}, {fine2, info2, "2 levels"}}
    [mesh_f, info, label] = pair{1}{:};
    u_fine = FEField(basis, L2proj2D(mesh_f, p, k));
    e = refinement_l2_error(u_coarse, coarse, u_fine, mesh_f, info);
    n_failed = n_failed + check(sprintf("exact polynomial, %s (%.2e)", label, e), ...
        e < 1e-12);
end

%% A field against itself at zero refinement ------------------------------
info0 = struct('parent_element', 1:size(coarse.elements, 2), ...
    'child_type', ones(1, size(coarse.elements, 2)), ...
    'A', repmat(eye(2), 1, 1, size(coarse.elements, 2)), ...
    'b', zeros(2, size(coarse.elements, 2)), ...
    'parent_face', 1:size(coarse.faces, 2), ...
    'alpha', ones(1, size(coarse.faces, 2)), ...
    'beta', zeros(1, size(coarse.faces, 2)), 'n_types', 1);
e_self = refinement_l2_error(u_coarse, coarse, u_coarse, coarse, info0);
n_failed = n_failed + check(sprintf("field against itself (%.2e)", e_self), ...
    e_self < 1e-13);

%% Known answer: a non-representable function -----------------------------
% sin is not in either space, so the two representations genuinely differ,
% and the distance must match what compute_l2_error reports for the same
% pair of fields evaluated on the fine mesh.
g = @(x, y) sin(3*x) .* cos(2*y);
v_coarse = FEField(basis, L2proj2D(coarse, g, k));
v_fine   = FEField(basis, L2proj2D(fine1, g, k));
e_cross  = refinement_l2_error(v_coarse, coarse, v_fine, fine1, info1);
% Triangle inequality against the two errors measured the ordinary way.
e_c = compute_l2_error(v_coarse, coarse, Datum(g));
e_f = compute_l2_error(v_fine, fine1, Datum(g));
n_failed = n_failed + check(sprintf( ...
    "cross error %.3e within [|%.3e-%.3e|, %.3e]", e_cross, e_c, e_f, e_c + e_f), ...
    e_cross <= e_c + e_f + 1e-12 && e_cross >= abs(e_c - e_f) - 1e-12);
n_failed = n_failed + check("cross error is not degenerate", e_cross > 1e-8);

%% Trace: exactness on the inherited skeleton -----------------------------
% A degree-k polynomial in the face parameter is in both trace spaces.
% The projected function must be the same in physical space, not the same
% function of the face parameter: a child face carries a different
% parametrization from its parent, so projecting the same q(t) on both
% would compare two different functions and prove nothing.
face_basis = LegendreBasis(k);
g_face = @(x, y) 0.7 - 1.3*x + 0.4*y + 0.2*x.*y;
uh_coarse = FEField(face_basis, face_dofs(coarse, face_basis, g_face, k));
uh_fine1  = FEField(face_basis, face_dofs(fine1,  face_basis, g_face, k));
e_tr = refinement_l2_error(uh_coarse, coarse, uh_fine1, fine1, info1);
n_failed = n_failed + check(sprintf("trace, exact polynomial (%.2e)", e_tr), ...
    e_tr < 1e-12);

%% The trace norm convention -----------------------------------------------
% At zero refinement parent and child are the same face, so the weighting
% (parent length as the norm weight, child length as the Jacobian) must
% collapse to exactly what FEField.compute_l2_error uses on the skeleton --
% the mesh-dependent ||h^(1/2) e||, i.e. its area.^2 scaling.
g2 = @(x, y) 0.1 + x - 0.3*y;
wh_coarse = FEField(face_basis, face_dofs(coarse, face_basis, g2, k));
info0_face = info0;
e_mine = refinement_l2_error(uh_coarse, coarse, wh_coarse, coarse, info0_face);
e_repo = compute_l2_error(uh_coarse - wh_coarse, coarse);
n_failed = n_failed + check(sprintf( ...
    "trace norm matches compute_l2_error (%.3e vs %.3e)", e_mine, e_repo), ...
    abs(e_mine - e_repo) <= 1e-12 * max(1, e_repo));

%% The element restriction ------------------------------------------------
% Passing every element must reproduce the unrestricted number exactly, and
% dropping elements must strictly reduce the (squared) integral. Guards the
% subset path used when an unfitted hierarchy is re-carved and the two
% meshes cover different domains.
e_all = refinement_l2_error(v_coarse, coarse, v_fine, fine1, info1, ...
    elements = 1:size(fine1.elements, 2));
n_failed = n_failed + check(sprintf("restriction to all == unrestricted (%.3e)", ...
    abs(e_all - e_cross)), abs(e_all - e_cross) < 1e-15 * max(1, e_cross));
half = 1:2:size(fine1.elements, 2);
e_half = refinement_l2_error(v_coarse, coarse, v_fine, fine1, info1, elements = half);
n_failed = n_failed + check(sprintf("restriction drops mass (%.3e < %.3e)", ...
    e_half, e_cross), e_half < e_cross && e_half > 0);

%% Does the self-error track the true error? -------------------------------
% The whole premise of a self-convergence study: ||u_j - u_L|| stands in for
% ||u_j - u|| while the reference is much finer. Tested where the true error
% is known -- L2 projection of a smooth function is the best approximation
% in each space, so both errors are computable and both must show k+1.
%
% This is what separates "the plumbing runs" from "the measure is right",
% and it also exhibits the reference bias: the level closest to the
% reference is the one whose self-error understates the truth.
g3 = @(x, y) exp(0.7*x) .* sin(2.3*y);
n_levels = 4;
hier = cell(1, n_levels);  hier{1} = StructuredSimplicialGrid(2, 2);
for j = 2:n_levels
    hier{j} = uniform_refine(hier{j-1});
end
proj = cellfun(@(m) FEField(basis, L2proj2D(m, g3, k)), hier, ...
    UniformOutput = false);

true_err = zeros(1, n_levels - 1);
self_err = zeros(1, n_levels - 1);
h_lvl    = zeros(1, n_levels - 1);
for j = 1:n_levels - 1
    [~, info_j] = uniform_refine(hier{j}, n_levels - j);
    true_err(j) = compute_l2_error(proj{j}, hier{j}, Datum(g3));
    self_err(j) = refinement_l2_error(proj{j}, hier{j}, ...
        proj{n_levels}, hier{n_levels}, info_j);
    h_lvl(j) = max(hier{j}.area);
end
rate = @(e) log(e(1:end-1) ./ e(2:end)) ./ log(h_lvl(1:end-1) ./ h_lvl(2:end));
r_true = rate(true_err);  r_self = rate(self_err);
fprintf("      true errors %s\n", sprintf("%.3e ", true_err));
fprintf("      self errors %s\n", sprintf("%.3e ", self_err));
fprintf("      true rates  %s\n", sprintf("%.2f ", r_true));
fprintf("      self rates  %s\n", sprintf("%.2f ", r_self));

n_failed = n_failed + check(sprintf("true error converges at k+1 (%.2f)", r_true(1)), ...
    abs(r_true(1) - (k + 1)) < 0.25);
n_failed = n_failed + check(sprintf("self rate matches true rate (%.2f vs %.2f)", ...
    r_self(1), r_true(1)), abs(r_self(1) - r_true(1)) < 0.25);
% The coarsest level is far from the reference, so the two errors should
% nearly coincide there.
n_failed = n_failed + check(sprintf("self error tracks true error (%.1f%%)", ...
    100*abs(self_err(1) - true_err(1))/true_err(1)), ...
    abs(self_err(1) - true_err(1)) < 0.05 * true_err(1));
% ... and the level adjacent to the reference must understate it, so a
% self-convergence study should exclude that level.
n_failed = n_failed + check("reference bias understates the finest level", ...
    self_err(end) < true_err(end));

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function dofs = face_dofs(mesh, basis, g, k)
% L2-project the physical function g onto the trace space of every face.
% On a straight face g is a polynomial of the same degree in the face
% parameter, so degree <= k is represented exactly on parent and child
% alike -- which is what makes the comparison return zero.
[t, w] = get_1d_reference_quadrature(k + 2, "Simplex");
t = reshape(t, 1, []);
w = reshape(w, [], 1);
n_faces = size(mesh.faces, 2);
xy = get_1d_physical(mesh, repmat(t(:), 1, n_faces), 1:n_faces);
G = reshape(g(xy(1, :, :), xy(2, :, :)), numel(t), n_faces);
P = basis.evaluate_basis(k, t(:));
dofs = (P' * (w .* P)) \ (P' * (w .* G));
end

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
