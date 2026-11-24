function e = refinement_l2_error(coarse_field, coarse_mesh, fine_field, fine_mesh, info, options)
%-------------------------------------------------------------------------
%Description: L2 distance between a discrete field and a discrete field on
%             a uniform refinement of its mesh (self-convergence studies).
%
%             The coarse field is evaluated at the fine mesh's quadrature
%             points, once per distinct child map (info.A / info.b from
%             uniform_refine), and the difference integrated there; exact
%             up to quadrature, with no interpolation or point location.
%
%             Volume and trace fields are both handled, dispatching on the
%             basis. The trace is measured on the inherited skeleton only.
%
%Usage:       [fine_mesh, info] = uniform_refine(coarse_mesh, 2);
%             e = refinement_l2_error(u_H, coarse_mesh, u_h, fine_mesh, info)
%
%See also: SimplicialGrid.uniform_refine, FEField.compute_l2_error
%-------------------------------------------------------------------------
arguments
    coarse_field FEField
    coarse_mesh  SimplicialGrid
    fine_field   FEField
    fine_mesh    SimplicialGrid
    info         struct
    % Quadrature rule index (exact to degree 2*order + 1; default k + 2).
    options.order double {mustBeScalarOrEmpty} = []
    % Restrict the integral to this subset of fine entities (elements, or
    % faces for a trace field). Empty = all of them. Use it when the fine
    % mesh covers more than the coarse one (e.g. a re-carved unfitted
    % hierarchy).
    options.elements double = []
end
order = options.order;
if isempty(order), order = coarse_field.basis.Order + 2; end

if coarse_field.basis.GeometricInputDimension == 1
    e = trace_error(coarse_field, coarse_mesh, fine_field, fine_mesh, info, ...
        order, options.elements);
else
    e = volume_error(coarse_field, coarse_mesh, fine_field, fine_mesh, info, ...
        order, options.elements);
end
end

%=========================================================================%
function e = volume_error(coarse_field, coarse_mesh, fine_field, fine_mesh, info, ...
    order, subset)
[rs, w] = get_2d_reference_quadrature(order, "Simplex");

n_fine = size(fine_mesh.elements, 2);
U_fine = flatten(fine_field.evaluate(rs, fine_mesh), numel(w), n_fine);

% volume is signed (shoelace)
measure = abs(fine_mesh.volume(:)).';
included = subset_mask(subset, n_fine);

e2 = 0;
for c = unique(info.child_type(included))
    j = find(info.child_type == c & included);
    K = info.parent_element(j);
    % same map for every element of this type: one broadcast call
    rs_coarse = info.A(:, :, j(1)) * rs + info.b(:, j(1));
    U_coarse = flatten(coarse_field.evaluate(rs_coarse, coarse_mesh), ...
        numel(w), size(coarse_mesh.elements, 2));

    D = U_fine(:, :, j) - U_coarse(:, :, K);
    per_element = reshape(sum(w .* D.^2, [1 2]), 1, []);
    e2 = e2 + sum(measure(j) .* per_element);
end
e = sqrt(e2);
end

%=========================================================================%
function e = trace_error(coarse_field, coarse_mesh, fine_field, fine_mesh, info, ...
    order, subset)
[t, w] = get_1d_reference_quadrature(order, "Simplex");
t = reshape(t, 1, []);
w = reshape(w, [], 1);

n_fine_faces = size(fine_mesh.faces, 2);
inherited = find(info.parent_face > 0 & subset_mask(subset, n_fine_faces));
U_fine = flatten(fine_field.evaluate(t, fine_mesh), numel(w), n_fine_faces);

% Mesh-dependent trace norm ||h^(1/2) e||: with unit-sum 1-D weights,
% sum_F |F|^2 w'e^2 = sum_F |F| int_F e^2. The weight h is the coarse face
% length; the Jacobian is the child face length.
parent_length = coarse_mesh.area(info.parent_face(inherited));
child_length  = fine_mesh.area(inherited);

% alpha, beta take four values (two halves x two stored directions)
maps = [info.alpha(inherited); info.beta(inherited)].';
[distinct, ~, which] = unique(maps, 'rows');

e2 = 0;
for m = 1:size(distinct, 1)
    sel = which == m;
    faces = inherited(sel);
    F = info.parent_face(faces);
    t_coarse = distinct(m, 1) * t + distinct(m, 2);
    U_coarse = flatten(coarse_field.evaluate(t_coarse, coarse_mesh), ...
        numel(w), size(coarse_mesh.faces, 2));

    D = U_fine(:, :, faces) - U_coarse(:, :, F);
    per_face = reshape(sum(w .* D.^2, [1 2]), 1, []);
    e2 = e2 + sum(parent_length(sel) .* child_length(sel) .* per_face);
end
e = sqrt(e2);
end

%=========================================================================%
function mask = subset_mask(subset, n_entities)
% Logical row over all fine entities; empty means "everything".
if isempty(subset)
    mask = true(1, n_entities);
else
    mask = false(1, n_entities);
    mask(subset) = true;
end
end

%=========================================================================%
function U = flatten(U, n_points, n_entities)
% FEField.evaluate returns (n_points x 1 x out1 x out2 x n_entities); the
% component structure is irrelevant to a Frobenius L2 norm, so collapse it.
U = reshape(U, n_points, [], n_entities);
end
