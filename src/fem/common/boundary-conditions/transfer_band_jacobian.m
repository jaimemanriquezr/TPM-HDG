function [J, info] = transfer_band_jacobian(mesh, mapping, options)
%TRANSFER_BAND_JACOBIAN  Signed area element of the transfer-path band.
%   J = TRANSFER_BAND_JACOBIAN(MESH, MAPPING, quadrature_face=s,
%   quadrature_path=t, tag=...) returns the signed area element of the band
%   swept by the transfer paths between MESH's tagged boundary and the true
%   boundary reached by MAPPING, per quadrature node, as an
%   (n_s x n_t x n_faces) array. Contracted with the tensorized face/path
%   weights it integrates over the band; contracted with 1 it gives the
%   signed area the computational domain is missing (positive) or covers in
%   excess (negative).
%
%   Both quadrature counts are required: a caller integrating a field over
%   the band must use the same rule it evaluates the field on.
%
%   The band is the ruled surface X(s,t) = (1-t) x(s) + t xbar(s), where x(s)
%   runs along the straight boundary face and xbar(s) = MAPPING(x(s)) along
%   the true boundary. With m = xbar - x and m' = xbar' - x',
%
%       X_s = x' + t m',   X_t = m,   det DX = x' x m + t (m' x m).
%
%   t = 0 is the mesh point and t = 1 the true boundary, as in the path
%   nodes of TransferEssentialBoundaryConditions. Since s, t run over [0,1]
%   and x' is the full chord vector, the tensorized weights contracted with
%   J give the signed band area directly.
%
%   Orientation: faces are traversed counterclockwise (chord = exterior
%   normal rotated +90 degrees), and J is negated, so a positive J means
%   area the computational domain is missing and a negative J area it
%   covers in excess. The band integral of 1 is then the signed area
%   difference area(true domain) - sum_K |K|, also where the discrete
%   boundary crosses the true one (folds of the ruled surface carry
%   opposite signs and cancel).
%
%   Name-value options:
%       tag               boundary tag bounding the band  (default "Dirichlet")
%       quadrature_face   s-nodes per face                (required)
%       quadrature_path   t-nodes per path                (required)
%       mapping_jacobian  optional handle returning D(MAPPING) as (2 x 2 x N)
%                         at (2 x N) points; when absent xbar' is obtained by
%                         central differences along the face
%       fd_step           step for those differences, in the face parameter
%                         s in [0,1]                      (default 1e-6)
%
%   The finite-difference xbar' limits the accuracy to about eps^(2/3);
%   supply mapping_jacobian when more is needed.
%
%   INFO reports .orientation (the factor applied, always -1), .net_gap
%   (the first-order signed gap sum(n.m), n the face exterior normal) and
%   .under_covers (true when net_gap > 0), plus .faces, .elements,
%   .sides, .m, the nodes and .weights (contract J with INFO.weights to
%   integrate over the band). INFO.x and INFO.xbar are the face and
%   true-boundary points of the oriented traversal; use them to build
%   X(s,t) = x + t*m rather than rebuilding from mesh.faces.
%
%   See also TRANSFERESSENTIALBOUNDARYCONDITIONS, IMMERSE_MESH.
arguments
    mesh                          SimplicialGrid
    mapping                 (1,1) {mustBeA(mapping, ["function_handle", "LevelSetCurve"])}
    options.tag             (1,1) string = "Dirichlet"
    options.quadrature_face (1,1) double          % required: see help
    options.quadrature_path (1,1) double          % required: see help
    options.mapping_jacobian      function_handle {mustBeScalarOrEmpty} = ...
                                  function_handle.empty
    options.fd_step         (1,1) double = 1e-6
end
if isa(mapping, "LevelSetCurve")
    mapping = @(X) mapping.project(X);
end

faces = get_face_by_tag(mesh, options.tag);
faces = faces(:).';
n_f = numel(faces);

[K, side] = get_element_by_face(mesh, faces);
if iscell(K)
    % a face on the tagged boundary borders exactly one element; the cell
    % form appears only for interior faces, so take the single entry
    K = cellfun(@(C) C(1), K);
    side = cellfun(@(C) C(1), side);
end

[s_nodes, w_s] = get_1d_reference_quadrature(options.quadrature_face);
[t_nodes, w_t] = get_1d_reference_quadrature(options.quadrature_path);
[s_nodes, w_s, t_nodes, w_t] = deal(s_nodes(:), w_s(:), t_nodes(:), w_t(:));
n_s = numel(s_nodes);
n_t = numel(t_nodes);

% ---- face geometry ---------------------------------------------------- %
normals = reshape(mesh.exterior_unit_normals(:, 1, :, :), 2, 3, []);
n_face = zeros(2, 1, n_f);
for b = 1:n_f
    n_face(:, 1, b) = normals(:, side(b), K(b));
end

X = reshape(mesh.coordinates(:, mesh.faces(1:2, :)), 2, 2, []);
X = X(:, :, faces);
chord = diff(X, 1, 2);                                   % (2 x 1 x n_f)

% mesh.faces stores vertices in global index order, so the raw chord
% direction is arbitrary per face. Traverse counterclockwise (chord = exterior
% normal rotated +90 degrees) so det DX does not alternate in sign.
flip = double(sum(chord .* [-n_face(2, :, :); n_face(1, :, :)], 1) < 0);
x_start = (1 - flip) .* X(:, 1, :) + flip .* X(:, 2, :);
chord = (1 - 2 * flip) .* chord;

% x(s) = x_start + s * chord, so chord is dx/ds on s in [0,1]
x_s = x_start + reshape(s_nodes, 1, n_s) .* chord;       % (2 x n_s x n_f)

xbar = apply_map(mapping, x_s);
m = xbar - x_s;

% ---- xbar' ------------------------------------------------------------ %
if isempty(options.mapping_jacobian)
    % central difference in s; a step of d in s is a step of d*chord in x.
    % Gauss nodes are strictly interior to [0,1], so s +/- d stays on the
    % face for any sane fd_step.
    d = options.fd_step;
    xbar_prime = (apply_map(mapping, x_s + d * chord) ...
                - apply_map(mapping, x_s - d * chord)) / (2 * d);
else
    D = options.mapping_jacobian(reshape(x_s, 2, []));    % (2 x 2 x n_s*n_f)
    chord_rep = reshape(repmat(chord, 1, n_s, 1), 2, 1, []);
    xbar_prime = reshape(pagemtimes(D, chord_rep), 2, n_s, n_f);
end
m_prime = xbar_prime - chord;

% ---- det DX = x' x m + t (m' x m) ------------------------------------- %
c0 = cross2(repmat(chord, 1, n_s, 1), m);                % (1 x n_s x n_f)
c1 = cross2(m_prime, m);                                 % (1 x n_s x n_f)

J = reshape(c0, n_s, 1, n_f) ...
  + reshape(c1, n_s, 1, n_f) .* reshape(t_nodes, 1, n_t);

% ---- orientation ------------------------------------------------------ %
% In the counterclockwise frame tangent x normal = -1, and m runs along +n
% when the paths leave the domain, so the raw determinant is negative where
% the domain is missing area. Negating makes a positive J mean missing area.
orientation = -1;
J = J * orientation;

signed_gap = sum(n_face .* m, 1);                        % (1 x n_s x n_f)
face_length = reshape(mesh.area(faces), 1, []);
net_gap = sum(reshape(signed_gap, n_s, n_f) .* w_s .* face_length, 'all');

info = struct( ...
    "orientation", orientation, "net_gap", net_gap, ...
    "under_covers", net_gap > 0, ...
    "faces", faces, "elements", K(:).', "sides", side(:).', ...
    "m", m, "x", x_s, "xbar", xbar, ...
    "s_nodes", s_nodes, "t_nodes", t_nodes, ...
    "s_weights", w_s, "t_weights", w_t, ...
    "weights", w_s * w_t.');   % (n_s x n_t), contract with J to integrate
end

% ----------------------------------------------------------------------- %
function Y = apply_map(f, X)
%APPLY_MAP  Evaluate a (2 x N) -> (2 x N) point map on a paged array.
sz = size(X);
Y = reshape(f(reshape(X, 2, [])), sz);
end

function c = cross2(a, b)
%CROSS2  Scalar cross product of paged 2-vectors, along dimension 1.
c = a(1, :, :) .* b(2, :, :) - a(2, :, :) .* b(1, :, :);
end
