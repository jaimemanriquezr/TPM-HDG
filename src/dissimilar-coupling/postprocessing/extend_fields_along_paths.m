function [x_band, u_band, sigma_band] = extend_fields_along_paths( ...
    problem, solution, mesh, partition, fractions, options)
%-------------------------------------------------------------------------
%Description: Extends the displacement (or velocity) and the stress off
%             the mesh, along the TPM transfer paths, to the points
%
%               x(s) = x + s m,   m = xbar - x,   s in [0, 1],
%
%             from the mesh-boundary point x (s = 0) to the interface
%             point xbar (s = 1). The fundamental theorem of calculus along
%             the path, with the path gradient reconstructed from the
%             mixed variables of the element K owning the face,
%
%               u(x(s)) = u_hat(x) + s int_0^1 grad u(x(s eta)) . m d eta,
%               grad u  = c1 sigma - c2 tr(sigma) I + gamma,
%               sigma(x(s)) = sigma_K(x(s)),
%
%             is the construction of transfer_trace_to_interface, which is
%             the s = 1 case. It serves both sides of the FSI coupling:
%
%               LinearElasticityProblem          c1 = 1/(2 mu),
%                                                c2 = lambda/(2 mu (2 lambda + 2 mu))
%               IncompressibleElasticityProblem  c1 = 1/(2 mu), c2 = 1/(4 mu)
%
%             (the fluid constants are those of FluidStructureCoupling.solve).
%-------------------------------------------------------------------------
arguments (Input)
    problem     (1,1)   FEProblem       % supplies the constitutive constants
    solution    (1,1)   FESolution      % post-solve fields (sigma, u, gamma, u_hat)
    mesh                SimplicialGrid  % mesh the solution lives on
    partition   (1,1)   struct          % this mesh's induced partition (get_partitions)
    fractions   (1,:)   double = linspace(0, 1, 5)   % the s values to sample
    options.Order double {mustBeScalarOrEmpty} = []
end
arguments (Output)
    x_band      (2,:,:,:)   double  % 2 x n_patch_points x n_fractions x n_patches
    u_band      (2,:,:,:)   double  % same layout
    sigma_band  (2,2,:,:,:) double  % 2 x 2 x n_patch_points x n_fractions x n_patches
end
k = options.Order;
if isempty(k), k = problem.Order; end

if isa(problem, "IncompressibleElasticityProblem")
    mu = problem.Data.viscosity;
    c1 = 1/(2*mu);
    c2 = 1/(4*mu);
else
    [lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
        problem.Data.young_modulus, problem.Data.poisson_ratio);
    c1 = 1/(2*mu);
    c2 = lambda/(2*mu*(2*lambda + 2*mu));
end

s_quad = size(partition.quadrature.local, 2);
n_frac = numel(fractions);
% grad u is a polynomial of degree k along a straight path: k + 1 Gauss
% points integrate it exactly; one more is margin.
[eta, w_eta] = lgwt(k + 2, 0, 1);
n_eta = numel(eta);
n_of_patches = size(partition.faces, 3);
G = get_affine_map(mesh, partition.elements, inverse = true);

sigma_dof = solution.Unknowns.field(1).dof;
gamma_dof = solution.Unknowns.field(3).dof;
if iscell(gamma_dof), gamma_dof = cat(1, gamma_dof{:}); end
u_hat = solution.Unknowns.hybrid(1).dof;
if ~iscell(u_hat)
    % flattened (2*d1, n_faces); restore the {u_x, u_y} view
    d1 = size(u_hat, 1) / 2;
    u_hat = {u_hat(1:d1, :), u_hat(d1 + 1:end, :)};
end

x_band = zeros(2, s_quad, n_frac, n_of_patches);
u_band = zeros(2, s_quad, n_frac, n_of_patches);
sigma_band = zeros(2, 2, s_quad, n_frac, n_of_patches);

for e = 1:n_of_patches
    K = partition.elements(e);
    F = partition.faces(e);
    x0 = partition.quadrature.physical(:, :, e);                 % 2 x s
    m = partition.quadrature.interface(:, :, e) - x0;            % 2 x s
    l = legendre_basis(k, partition.quadrature.local(1, :, e));  % s x d1
    u0 = l * [u_hat{1}(:, F), u_hat{2}(:, F)];                   % s x 2

    for j = 1:n_frac
        s = fractions(j);
        xs = x0 + s*m;
        sigma_s = local_fields(mesh, K, G(:, :, e), xs, k, sigma_dof, gamma_dof);

        % Gauss points along [0, s]: 2 x s x n_eta, flattened point-first
        xq = reshape(x0, 2, s_quad) + s*reshape(m, 2, s_quad) .* reshape(eta, 1, 1, n_eta);
        [sigma_q, gamma_q] = local_fields(mesh, K, G(:, :, e), reshape(xq, 2, []), ...
            k, sigma_dof, gamma_dof);
        trace_q = sigma_q(:, 1, 1) + sigma_q(:, 2, 2);
        grad_u = c1*sigma_q + gamma_q;                           % N x 2 x 2
        grad_u(:, 1, 1) = grad_u(:, 1, 1) - c2*trace_q;
        grad_u(:, 2, 2) = grad_u(:, 2, 2) - c2*trace_q;
        mq = reshape(repmat(m, [1, 1, n_eta]), 2, []).';         % N x 2
        grad_u_m = [grad_u(:, 1, 1).*mq(:, 1) + grad_u(:, 1, 2).*mq(:, 2), ...
                    grad_u(:, 2, 1).*mq(:, 1) + grad_u(:, 2, 2).*mq(:, 2)];
        grad_u_m = reshape(grad_u_m, s_quad, n_eta, 2);
        path_integral = s * reshape(sum(reshape(w_eta, 1, []) .* grad_u_m, 2), s_quad, 2);

        x_band(:, :, j, e) = xs;
        u_band(:, :, j, e) = (u0 + path_integral).';
        sigma_band(:, :, :, j, e) = permute(sigma_s, [2 3 1]);
    end
end
end

%=========================================================================%
function [sigma, gamma] = local_fields(mesh, K, G_e, X, k, sigma_dof, gamma_dof)
% sigma and gamma of element K at the physical points X (2 x N), which may
% lie outside K: the element polynomials extrapolate. N x 2 x 2 each.
reference = reshape(mesh.get_2d_local(reshape(X, 2, [], 1), K), 2, []);
P = PEERS_basis(k, reference(1, :).', reference(2, :).', G_e);      % N x nb x 2 x 2
sigma = reshape(tensorprod(P, sigma_dof(:, K), 2, 1), [], 2, 2);
if nargout > 1
    A = asym_tensor_dubiner_basis(k, reference(1, :).', reference(2, :).');
    gamma = reshape(tensorprod(A, gamma_dof(:, K), 2, 1), [], 2, 2);
end
end
