function [x_interface, t_interface] = transfer_traction_to_interface( ...
    solution, mesh, partition, options)
%-------------------------------------------------------------------------
%Description: Evaluates the transferred traction E(sigma)n on the interface
%             Sigma from the recovered volume stress, by the same
%             extrapolated pointwise evaluation the coupled solver uses in
%             its T1 (traction transmission) rows.
%
%             Unlike the displacement transfer (transfer_trace_to_interface,
%             a path integral of the reconstructed gradient), the traction
%             transfer is a pointwise extrapolation: the PEERS stress basis
%             of the patch's element is evaluated directly at the interface
%             quadrature points (reference coordinates
%             partition.quadrature.interface_2d, which lie outside the
%             element wherever there is a gap) and contracted with the
%             supplied normal.
%
%             The caller chooses the normal convention. The solver's T1 row
%             uses n_F = -partition.quadrature.normals on the fluid side
%             and n_S = +partition.quadrature.normals on the solid side.
%-------------------------------------------------------------------------
arguments (Input)
    solution    (1,1)   FESolution     % post-solve fields (sigma in volume(1))
    mesh                SimplicialGrid  % mesh the solution lives on
    partition   (1,1)   struct          % this mesh's induced partition (get_partitions)
    options.Order (1,1) double
    % Signed normals at the interface quadrature points, 2 x s_quad x
    % n_patches. Defaults to +partition.quadrature.normals (the canonical
    % interface normal; negate for the fluid side of the T1 convention).
    options.normals double = []
end
arguments (Output)
    x_interface (2,:,:) double  % interface points, 2 x s_quad x n_patches
    t_interface (2,:,:) double  % traction E(sigma)n there, same layout
end
k = options.Order;
normals = options.normals;
if isempty(normals), normals = partition.quadrature.normals; end

s_quad = size(partition.quadrature.local, 2);
n_of_patches = size(partition.faces, 3);
G = get_affine_map(mesh, partition.elements, inverse = true);

% sigma: Dubiner component blocks + curl-curl bubble block, stacked in the
% column order PEERS_basis produces (same recovery as the displacement
% transfer).
sigma_dof = solution.Unknowns.field(1).dof;

[x_interface, t_interface] = deal(zeros(2, s_quad, n_of_patches));

for e = 1:n_of_patches
    K = partition.elements(e);

    s_iface = partition.quadrature.interface_2d(1, :, e);
    t_iface = partition.quadrature.interface_2d(2, :, e);

    P = PEERS_basis(k, s_iface(:), t_iface(:), G(:, :, e));   % s x nb x 2 x 2
    sigma = tensorprod(P, sigma_dof(:, K), 2, 1);             % s x 2 x 2

    n = normals(:, :, e);                                     % 2 x s
    t = zeros(2, s_quad);
    t(1, :) = squeeze(sigma(:, 1, 1)).' .* n(1, :) + squeeze(sigma(:, 1, 2)).' .* n(2, :);
    t(2, :) = squeeze(sigma(:, 2, 1)).' .* n(1, :) + squeeze(sigma(:, 2, 2)).' .* n(2, :);

    t_interface(:, :, e) = t;
    x_interface(:, :, e) = partition.quadrature.interface(:, :, e);
end
end
