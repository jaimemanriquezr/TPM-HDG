function J = get_sliver_jacobian(interface, mesh, partition, options)
%GET_SLIVER_JACOBIAN Signed area element of the TPM gap sliver.
%   J = GET_SLIVER_JACOBIAN(interface, mesh, partition, patch_weight=w,
%   quadrature_patch=s, quadrature_path=t) returns the signed area element
%   of the sliver between the mesh's discrete interface and Sigma, per
%   patch quadrature node, as an s x t x n array (patch node,
%   transfer-path node, patch). All three name-value arguments are
%   required. Multiplied by the tensorized patch/path
%   quadrature weights it integrates over the sliver; summed it gives the
%   signed area the mesh's computational domain is missing (positive,
%   paths extrapolate out of the mesh domain) or covers in excess
%   (negative, paths run into it, e.g. the outer mesh of a fitted pair
%   around a convex Sigma, whose chords cut into the inner subdomain), on
%   either side of Sigma and for either patch weight.
%
%   The sliver is the ruled surface X(s,t) = (1-t) xbar(s) + t x(s) with
%   det DX = -( x' x m + (1-t) m' x m ), m = xbar - x. Each patch is
%   traversed counterclockwise around the mesh's own domain (chord = the
%   face's exterior normal rotated +90 degrees), so x' x n < 0 and the
%   determinant is positive where m runs along the exterior normal. The
%   determinant stays signed pointwise: where the discrete boundary weaves
%   across Sigma the folds cancel exactly, while a per-node sign
%   double-counts the crossed lobes.
%
%   J = GET_SLIVER_JACOBIAN(..., patch_weight="rectangular") returns the
%   theory-conformal rectangular weight |patch| * |m(s)| * sign(n . m),
%   n the face's exterior unit normal, instead: valid when the transfer
%   paths stay close to the face normals (theta -> 0 with h), where both
%   weights are norm-equivalent. With oblique paths (e.g. the canonical
%   vertical pairing, theta = O(1)) it carries the O(delta) alpha-row
%   defect.

arguments
    interface AbstractCurve  % interface Sigma (provides tangent)
    mesh      SimplicialGrid         % mesh whose boundary bounds the sliver
    partition struct                 % induced partition from get_partitions
    % No defaults: patch_weight is sign-critical and the quadrature orders
    % must match the solver's. Checked explicitly below.
    options.patch_weight     (1,1) string = missing
    options.quadrature_patch (1,1) double = nan  % s-nodes per patch (match solver)
    options.quadrature_path  (1,1) double = nan  % t-nodes per path (match solver)
end
assert(~ismissing(options.patch_weight), "get_sliver_jacobian:missingArgument", ...
    "patch_weight is required (sign-critical): pass ""jacobian"" or ""rectangular"".")
mustBePatchWeight(options.patch_weight)
assert(~isnan(options.quadrature_patch) && ~isnan(options.quadrature_path), ...
    "get_sliver_jacobian:missingArgument", ...
    "quadrature_patch and quadrature_path are required and must match the solver's orders.")
s_quad  = options.quadrature_patch;
tp_quad = options.quadrature_path;
xiS = lgwt(s_quad, 0, 1);
tTP = lgwt(tp_quad, 0, 1);

m = partition.quadrature.interface - partition.quadrature.physical;   % 2 x s x n

% Exterior unit normal of the mesh face carrying each patch
normals = reshape(mesh.exterior_unit_normals(:, 1, :, :), 2, []);
n_face = reshape(normals(:, 3*(partition.elements(:) - 1) + partition.side_number(:)), 2, 1, []);

if options.patch_weight == "rectangular"
    % The rectangular weight carries no signed determinant, so the
    % per-node sign is its only signed content and stays pointwise:
    % positive where the path leaves the mesh domain.
    rect_core = reshape(partition.area, 1, 1, []) .* vecnorm(m, 2, 1) ...
        .* sign(sum(n_face .* m, 1));
    J = repmat(permute(rect_core, [2 1 3]), 1, tp_quad, 1);
    return
end

s_ep = partition.vertices_interface;            % 1 x 2 x n, global params
s_glob = s_ep(1, 1, :) + reshape(xiS, 1, [], 1) .* (s_ep(1, 2, :) - s_ep(1, 1, :));
xbar_prime = reshape(interface.tangent(reshape(s_glob, 1, [])), 2, s_quad, []) ...
    .* (s_ep(1, 2, :) - s_ep(1, 1, :));          % 2 x s x n
X_ends = get_1d_physical(mesh, partition.vertices, partition.faces);
x_prime = X_ends(:, 2, :) - X_ends(:, 1, :);    % 2 x 1 x n (patch chord)
% Counterclockwise traversal (see help above). Reversing s -> 1-s negates
% x' and xbar' at the same physical point and leaves m unchanged, so the
% quadrature nodes keep their order.
flip = sum(x_prime .* [-n_face(2, :, :); n_face(1, :, :)], 1) < 0;
x_prime = (1 - 2 * flip) .* x_prime;
xbar_prime = (1 - 2 * flip) .* xbar_prime;
m_prime = xbar_prime - x_prime;
cross2 = @(a, b) a(1,:,:).*b(2,:,:) - a(2,:,:).*b(1,:,:);
c0_sliver = cross2(repmat(x_prime, 1, s_quad, 1), m);
c1_sliver = cross2(m_prime, m);
J = -(reshape(c0_sliver, s_quad, 1, []) + ...
    reshape(c1_sliver, s_quad, 1, []) .* reshape(1 - tTP, 1, tp_quad, 1));   % s x t x n
end
