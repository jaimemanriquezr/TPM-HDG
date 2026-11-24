function path_integral = path_gradient_integral(k, partition, e, G_e, sigma_dof, gamma_dof, c1, c2)
%-------------------------------------------------------------------------
%Description: Path integral of the reconstructed gradient on patch e,
%
%               int_0^1 grad u(x(xi)) . m dxi,   m = xbar - x,
%
%             at the patch quadrature points, with the path gradient
%             reconstructed from the extrapolated mixed variables of the
%             element owning the patch through the PEERS constitutive law
%
%               grad u = c1 sigma - c2 tr(sigma) I + gamma.
%
%             Quadrature sizes are read off the partition so that the
%             transfer uses exactly the rule the partition was built with.
%
%See also: transfer_trace_to_interface, trace_extension_coefficients
%-------------------------------------------------------------------------
arguments (Input)
    k           (1,1)   double
    partition   (1,1)   struct  % induced partition (get_partitions)
    e           (1,1)   double  % patch index
    G_e         (2,2)   double  % inverse affine map of the patch's element
    sigma_dof   (:,1)   double  % PEERS stress dofs of that element
    gamma_dof   (:,1)   double  % rotation dofs of that element
    c1          (1,1)   double
    c2          (1,1)   double
end
arguments (Output)
    path_integral (:,2) double  % s x 2
end
s_quad  = size(partition.quadrature.local, 2);
tp_quad = size(partition.quadrature.path_2d, 4);
[~, wTP] = lgwt(tp_quad, 0, 1);

% Reference coordinates of the path quadrature points, laid out
% (tp x s) so that the linear index matches the basis evaluation.
s_conn = reshape(partition.quadrature.path_2d(1, :, e, :), s_quad, tp_quad).';
t_conn = reshape(partition.quadrature.path_2d(2, :, e, :), s_quad, tp_quad).';

P = PEERS_basis(k, s_conn(:), t_conn(:), G_e);
P = reshape(P, [tp_quad, s_quad, size(P, 2:4)]);            % tp x s x nb x 2 x 2
A = asym_tensor_dubiner_basis(k, s_conn(:), t_conn(:));
A = reshape(A, [tp_quad, s_quad, size(A, 2:4)]);            % tp x s x d2 x 2 x 2

sigma = tensorprod(P, sigma_dof, 3, 1);                     % tp x s x 2 x 2
gamma = tensorprod(A, gamma_dof, 3, 1);                     % tp x s x 2 x 2

% grad u = c1 sigma - c2 tr(sigma) I + gamma
trace_sigma = sigma(:, :, 1, 1) + sigma(:, :, 2, 2);
grad_u = c1*sigma + gamma;
grad_u(:, :, 1, 1) = grad_u(:, :, 1, 1) - c2*trace_sigma;
grad_u(:, :, 2, 2) = grad_u(:, :, 2, 2) - c2*trace_sigma;

% Path vector m = xbar - x, unnormalized (carries the path length).
m = partition.quadrature.interface(:, :, e) - partition.quadrature.physical(:, :, e);
% m.' first: the component index must land in dim 4 (the contracted
% tensor index), not be folded in column-major with the patch index.
grad_u_m = sum(grad_u .* reshape(m.', [1, s_quad, 1, 2]), 4);  % tp x s x 2
path_integral = reshape(sum(wTP .* grad_u_m, 1), s_quad, 2);
end
