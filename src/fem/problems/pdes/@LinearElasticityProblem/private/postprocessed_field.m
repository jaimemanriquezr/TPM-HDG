function U_star = postprocessed_field(obj, mesh, k, U, sigma, gamma)
% Element-wise P_k+1 reconstruction of a displacement-type field.
%
%   U_star = postprocessed_field(obj, mesh, k, U, sigma, gamma)
%   U_star = postprocessed_field(obj, mesh, k, U, L, [])
%
% U (2*d2 x N, Dubiner of degree k) with its PEERS stress and rotation
% dofs; U_star (2*d2' x N, Dubiner of degree k + 1, d2' = (k+2)(k+3)/2)
% solves on each element, per component c,
%   (grad u*_c, grad w)_K = (q_c, grad w)_K   for all w in P_k+1,
%   mean_K u*_c = mean_K u_c,
% with q = c1 sigma - c2 tr(sigma) I + gamma the reconstructed gradient
% (row c is grad u_c), the local postprocessing of the HDG displacement.
% With gamma = [], the fifth argument holds the PEERS dofs of the gradient
% q itself (e.g. the discrete gradient L_h of a hyperelastic solid).
% The Dubiner basis is orthogonal and hierarchical, so the mean is the
% first coefficient and the P_k projection of u* is its first d2 modes.
N = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
d2p = nchoosek(k + 3, 2);
[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    obj.Data.young_modulus, obj.Data.poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(2*mu*(2*lambda + 2*mu));

[rs, w] = get_2d_reference_quadrature(k + 1);        % exact to 2k + 3
nq = numel(w);
G = get_affine_map(mesh, inverse = true);
[~, dPr, dPs] = dubiner_basis(k + 1, rs);           % nq x d2p
g = @(i, j) reshape(G(i, j, :), 1, 1, N);
dPx = dPr.*g(1, 1) + dPs.*g(2, 1);                  % nq x d2p x N
dPy = dPr.*g(1, 2) + dPs.*g(2, 2);

% q at the quadrature points: nq x 2 x 2 x N
P = PEERS_basis(k, rs(1, :), rs(2, :), G);          % nq x nb x 2 x 2 x N
nb = size(P, 2);
P = reshape(permute(P, [1 3 4 2 5]), 4*nq, nb, N);
s = reshape(pagemtimes(P, reshape(sigma, nb, 1, N)), nq, 2, 2, N);
if isempty(gamma)
    q = s;
else
    A = asym_tensor_dubiner_basis(k, rs(1, :), rs(2, :));   % nq x d2 x 2 x 2
    A = reshape(permute(A, [1 3 4 2]), 4*nq, d2);
    r = reshape(A*gamma, nq, 2, 2, N);
    q = c1*s + r;
    trace_s = s(:, 1, 1, :) + s(:, 2, 2, :);
    q(:, 1, 1, :) = q(:, 1, 1, :) - c2*trace_s;
    q(:, 2, 2, :) = q(:, 2, 2, :) - c2*trace_s;
end

stiffness = pagemtimes(w.*dPx, 'transpose', dPx, 'none') ...
          + pagemtimes(w.*dPy, 'transpose', dPy, 'none');    % d2p x d2p x N
free = 2:d2p;
U_star = zeros(2*d2p, N);
for c = 1:2
    load = pagemtimes(w.*dPx, 'transpose', reshape(q(:, c, 1, :), nq, 1, N), 'none') ...
         + pagemtimes(w.*dPy, 'transpose', reshape(q(:, c, 2, :), nq, 1, N), 'none');
    modes = pagemldivide(stiffness(free, free, :), load(free, :, :));
    rows = (c - 1)*d2p;
    U_star(rows + 1, :) = U((c - 1)*d2 + 1, :);
    U_star(rows + free, :) = reshape(modes, d2p - 1, N);
end
end
