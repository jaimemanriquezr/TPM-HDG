function U_hat = shifted_skeleton(obj, mesh, k, U_hat, sigma, gamma, shift)
% Face-wise L2 projection of translated skeleton traces.
%
%   U_hat = shifted_skeleton(obj, mesh, k, U_hat, sigma, gamma, shift)
%
% For each field i (cells U_hat{i}, 2*d1 x F, with its PEERS stress and
% rotation dofs sigma{i}, gamma{i}), the trace at x - shift is reached
% from the face point x along the segment to x - shift with the
% reconstructed gradient grad u = c1 sigma - c2 tr(sigma) I + gamma of
% one element K owning the face (K's polynomials, extrapolated beyond K),
%   u_hat(x - shift) ~ u_hat(x) - int_0^1 grad u(x - s shift) shift ds,
% as the TPM transfer extends a trace along its path, and projected onto
% the face's Legendre basis.
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
N = size(mesh.elements, 2);
F = size(mesh.faces, 2);
n_fields = numel(U_hat);
[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    obj.Data.young_modulus, obj.Data.poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(2*mu*(2*lambda + 2*mu));

[t, w] = get_1d_reference_quadrature(k + 1);        % k + 2 points on [0, 1]
nt = numel(t);
[s, ws] = lgwt(k + 2, 0, 1);
ns = numel(s);
l = legendre_basis(k, t);                           % nt x d1
projector = (l.'*(w.*l)) \ (w.*l).';               % d1 x nt

% each face is owned by one of its elements
owner = zeros(F, 1);
owner(mesh.boundaries(:)) = repelem(1:N, 3);
owned = accumarray(owner, (1:F).', [N, 1], @(v) {v});
G = get_affine_map(mesh, inverse = true);
[~, x3] = get_affine_map(mesh);
x_1 = mesh.coordinates(1:2, mesh.faces(1, :));
x_2 = mesh.coordinates(1:2, mesh.faces(2, :));

correction = zeros(2*d1, F, n_fields);
for K = find(~cellfun(@isempty, owned)).'
    faces = owned{K}.';
    % points x - s*shift: nt x ns x numel(faces)
    x = reshape(x_1(:, faces), 2, 1, 1, []) + reshape(x_2(:, faces) - x_1(:, faces), 2, 1, 1, []) ...
        .* reshape(t, 1, nt) - shift .* reshape(s, 1, 1, ns);
    r = G(:, :, K)*(reshape(x, 2, []) - x3(:, K));
    P = PEERS_basis(k, r(1, :), r(2, :), G(:, :, K));           % n x nb x 2 x 2
    A = asym_tensor_dubiner_basis(k, r(1, :), r(2, :));          % n x d2 x 2 x 2
    for i = 1:n_fields
        sigma_x = tensorprod(P, sigma{i}(:, K), 2, 1);           % n x 2 x 2
        gamma_x = tensorprod(A, gamma{i}(1:d2, K), 2, 1);
        grad_u = c1*sigma_x + gamma_x;
        trace_sigma = sigma_x(:, 1, 1) + sigma_x(:, 2, 2);
        grad_u(:, 1, 1) = grad_u(:, 1, 1) - c2*trace_sigma;
        grad_u(:, 2, 2) = grad_u(:, 2, 2) - c2*trace_sigma;
        % -grad u shift, integrated over s: nt x numel(faces) x 2
        grad_u_shift = -(grad_u(:, :, 1)*shift(1) + grad_u(:, :, 2)*shift(2));
        path = reshape(grad_u_shift, nt, ns, numel(faces), 2);
        path = reshape(sum(path .* reshape(ws, 1, ns), 2), nt, numel(faces), 2);
        correction(:, faces, i) = [projector*path(:, :, 1); projector*path(:, :, 2)];
    end
end
for i = 1:n_fields
    U_hat{i} = U_hat{i} + correction(:, :, i);
end
end
