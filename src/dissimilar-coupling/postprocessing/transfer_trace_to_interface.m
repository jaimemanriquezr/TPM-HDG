function [x_interface, u_interface] = transfer_trace_to_interface( ...
    problem, solution, mesh, partition, options)
%-------------------------------------------------------------------------
%Description: Evaluates the elasticity displacement on the interface Sigma
%             by transferring the computed trace off the mesh boundary
%             along the TPM transfer paths.
%
%             With a gap, the solid mesh boundary sits at a distance from
%             Sigma, so reading the skeleton trace u_hat on the interface
%             faces returns u at the wrong points (an O(gap) error that
%             does not vanish under h-refinement). The transfer path
%             x(xi) = xbar + xi (x - xbar) from the interface point xbar to
%             the mesh point x carries the fundamental theorem of calculus
%
%               u(xbar) = u_hat(x) + int_0^1 grad u(x(xi)) . m dxi,
%               m = xbar - x,
%
%             and the path gradient is reconstructed from the extrapolated
%             mixed variables through the constitutive law of the PEERS
%             formulation (zeta^-1 sigma - grad u + gamma = 0):
%
%               grad u = c1 sigma - c2 tr(sigma) I + gamma,
%               c1 = 1/(2 mu),  c2 = lambda/(2 mu (2 lambda + 2 mu)).
%
%             The solver applies the same reconstruction in its T2 row;
%             here it is contracted with the recovered local dofs instead
%             of the local solvers.
%-------------------------------------------------------------------------
arguments (Input)
    problem     (1,1)   LinearElasticityProblem % supplies the Lame parameters
    solution    (1,1)   FESolution             % post-solve fields (sigma, u, gamma, u_hat)
    mesh                SimplicialGrid          % mesh the solution lives on
    partition   (1,1)   struct                  % this mesh's induced partition (get_partitions)
    options.Order double {mustBeScalarOrEmpty} = []
end
arguments (Output)
    x_interface (2,:)   double  % interface points, 2 x (n_patch_points * n_patches)
    u_interface (2,:)   double  % displacement there, same column ordering
end
k = options.Order;
if isempty(k), k = problem.Order; end

[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    problem.Data.young_modulus, problem.Data.poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(2*mu*(2*lambda + 2*mu));

s_quad = size(partition.quadrature.local, 2);

n_of_patches = size(partition.faces, 3);
G = get_affine_map(mesh, partition.elements, inverse = true);

% sigma: 4 Dubiner component blocks + (k >= 1) the curl-curl bubble block,
% stacked in the column order PEERS_basis produces. gamma: one scalar
% Dubiner block, expanded to a 2x2 skew tensor by the basis itself.
sigma_dof = solution.Unknowns.field(1).dof;
gamma_dof = solution.Unknowns.field(3).dof;
if iscell(gamma_dof), gamma_dof = cat(1, gamma_dof{:}); end
u_hat = solution.Unknowns.hybrid(1).dof;
if ~iscell(u_hat)
    % assign_unknowns flattens the per-component cell into a
    % (2*d1, n_faces) matrix (components stacked per face column);
    % restore the {u_x, u_y} view this function reads.
    d1 = size(u_hat, 1) / 2;
    u_hat = {u_hat(1:d1, :), u_hat(d1 + 1:end, :)};
end

[x_interface, u_interface] = deal(zeros(2, s_quad, n_of_patches));

for e = 1:n_of_patches
    K = partition.elements(e);
    F = partition.faces(e);

    % int_0^1 grad u(x(xi)) . m dxi, grad u from (sigma, gamma) of K
    path_integral = path_gradient_integral(k, partition, e, G(:, :, e), ...
        sigma_dof(:, K), gamma_dof(:, K), c1, c2);                % s x 2

    % Skeleton trace at the same patch parameters, on the interface face.
    l = legendre_basis(k, partition.quadrature.local(1, :, e));  % s x d1
    u_boundary = l * [u_hat{1}(:, F), u_hat{2}(:, F)];           % s x 2

    u_interface(:, :, e) = (u_boundary + path_integral).';
    x_interface(:, :, e) = partition.quadrature.interface(:, :, e);
end

x_interface = reshape(x_interface, 2, []);
u_interface = reshape(u_interface, 2, []);
end
