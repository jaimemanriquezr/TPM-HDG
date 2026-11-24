function [u_star, divergence] = postprocess_velocity(mesh, k, u_hat, u_h, L_h)
arguments
    mesh SimplicialGrid
    k (1,1) double {mustBePositive}   % degree of the HDG solution (k >= 1)
    u_hat (2,1) cell   % {ux_hat; uy_hat}: velocity trace, d1 x Nfaces each (Legendre)
    u_h (2,1) cell     % {u1; u2}: velocity, d2 x Nelts each (Dubiner)
    L_h (2,2) cell     % {L11, L12; L21, L22}: velocity gradient, d2 x Nelts each
end
%-------------------------------------------------------------------------
%Description: Element-by-element divergence-free postprocessing of the HDG
%             velocity (Cockburn, Gopalakrishnan, Nguyen, Peraire, Sayas
%             2011, eq. 2.9, two-dimensional form): u* in P_{k+1}(K)^2 with
%               <(u* - uhat).n, mu>_F = 0                    mu in P_k(F)
%               <d_t(u*.n) - t.({L}^T n), d_t mu>_F = 0      mu in P_{k+1}(F)^perp
%               (u* - u_h, grad w)_K = 0                     w in P_k(K)
%               (curl u* - (L21 - L12), w b_K)_K = 0         w in P_{k-1}(K)
%             where {L} is the face average of L_h from both sides, t the
%             unit tangent, b_K the cubic bubble. u* is H(div)-conforming,
%             exactly divergence-free and converges at order k+2 for k >= 1.
%             L_h is the velocity gradient of the gradient-based solver or
%             its reconstruction from a stress-based one (see
%             gradient_from_stress).
%Output:    - u_star: {u1*; u2*}, dofs in dubiner_basis(k+1), n x Nelts
%           - divergence: max over elements of ||div u*||_{L2(K)}
%-------------------------------------------------------------------------
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
dm = nchoosek(k + 1, 2);          % dim P_{k-1}
n = nchoosek(k + 3, 2);           % dim P_{k+1}
rows_a = 3*d1;  rows_b = 3;  rows_c = d2 - 1;  rows_d = dm;
assert(rows_a + rows_b + rows_c + rows_d == 2*n, "postprocess_velocity:count");

A = zeros(2*n, 2*n, N_ELEMENTS);
rhs = zeros(2*n, 1, N_ELEMENTS);
[F, xy_0] = get_affine_map(mesh, output="array");
det_F = reshape(F(1,1,:).*F(2,2,:) - F(1,2,:).*F(2,1,:), 1, 1, []);
grad_phys = @(pr, ps) deal( ...
    (reshape(F(2,2,:), 1, 1, []) .* pr - reshape(F(2,1,:), 1, 1, []) .* ps) ./ det_F, ...
    (-reshape(F(1,2,:), 1, 1, []) .* pr + reshape(F(1,1,:), 1, 1, []) .* ps) ./ det_F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VOLUME ROWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, rs, w] = quadpoints_elements_2d(mesh, k + 1, output="cell");
[P1, P1r, P1s] = dubiner_basis(k + 1, rs{:});
[P1x, P1y] = grad_phys(P1r, P1s);                       % nq x n x N
[Pk, Pkr, Pks] = dubiner_basis(k, rs{:});
[Pkx, Pky] = grad_phys(Pkr, Pks);                       % nq x d2 x N
vol = reshape(mesh.volume, 1, 1, []);
u1 = Pk*u_h{1};  u2 = Pk*u_h{2};                        % nq x N
u1 = reshape(u1, [], 1, N_ELEMENTS);  u2 = reshape(u2, [], 1, N_ELEMENTS);

% (u* - u_h, grad w) = 0, w in P_k without the constant
rc = rows_a + rows_b + (1:rows_c);
A(rc, 1:n, :)     = vol .* pagemtimes(w .* Pkx(:, 2:end, :), 'transpose', P1, 'none');
A(rc, n+1:2*n, :) = vol .* pagemtimes(w .* Pky(:, 2:end, :), 'transpose', P1, 'none');
rhs(rc, 1, :) = vol .* (pagemtimes(w .* Pkx(:, 2:end, :), 'transpose', u1, 'none') ...
                      + pagemtimes(w .* Pky(:, 2:end, :), 'transpose', u2, 'none'));

% (curl u* - (L21 - L12), w b_K) = 0, w in P_{k-1}, curl u = d_x u2 - d_y u1
rd = rows_a + rows_b + rows_c + (1:rows_d);
Pm = dubiner_basis(k - 1, rs{:});
bK = ((1 - rs{1} - rs{2}) .* rs{1} .* rs{2}).';        % nq x 1
wb = w .* bK .* Pm;                                     % nq x dm
omega = reshape(Pk*(L_h{2,1} - L_h{1,2}), [], 1, N_ELEMENTS);
A(rd, 1:n, :)     = -vol .* pagemtimes(wb, 'transpose', P1y, 'none');
A(rd, n+1:2*n, :) =  vol .* pagemtimes(wb, 'transpose', P1x, 'none');
rhs(rd, 1, :) = vol .* pagemtimes(wb, 'transpose', omega, 'none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FACE ROWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s, w1] = get_1d_reference_quadrature(k + 1);
nq = length(w1);
q2d = [[0*s(:), 1 - s(:)]; [1 - s(:), 0*s(:)]; [1 - s(:), s(:)]];
[P1e, P1er, P1es] = dubiner_basis(k + 1, q2d(:,1), q2d(:,2));
[P1ex, P1ey] = grad_phys(P1er, P1es);                   % 3nq x n x N
Pke = dubiner_basis(k, q2d(:,1), q2d(:,2));             % 3nq x d2
[l, dl] = legendre_basis(k + 1, [s(:); 1 - s(:)]);      % 2nq x (k+2)
xy_e = pagemtimes(F, q2d.') + reshape(xy_0, 2, 1, []);  % 2 x 3nq x N

area = mesh.area(mesh.boundaries);                      % 3 x N
oriented_normals = mesh.normals .* kron(mesh.orientation, [1;1]);
for e = 1:3
    nodes_e = (e-1)*nq + (1:nq);
    f = mesh.boundaries(e, :);
    len = reshape(area(e, :), 1, 1, []);
    nx = reshape(oriented_normals(2*e-1, :) ./ area(e, :), 1, 1, []);
    ny = reshape(oriented_normals(2*e, :) ./ area(e, :), 1, 1, []);
    tx = -ny;  ty = nx;

    % face basis with the element's orientation of the face parameter
    l_e = zeros(nq, k + 2, N_ELEMENTS);  dl_e = l_e;
    for mu = 1:2
        is_mu = reshape(mesh.permutation(e, :) == mu, 1, 1, []);
        l_e  = l_e  + is_mu .* l((mu-1)*nq + (1:nq), :);
        dl_e = dl_e + is_mu .* dl((mu-1)*nq + (1:nq), :);
    end

    % <(u* - uhat).n, mu> = 0, mu in P_k(F)
    ra = (e-1)*d1 + (1:d1);
    wl = w1 .* l_e(:, 1:d1, :);                         % nq x d1 x N
    A(ra, 1:n, :)     = len .* nx .* pagemtimes(wl, 'transpose', P1e(nodes_e, :), 'none');
    A(ra, n+1:2*n, :) = len .* ny .* pagemtimes(wl, 'transpose', P1e(nodes_e, :), 'none');
    ux = pagemtimes(l_e(:, 1:d1, :), reshape(u_hat{1}(:, f), d1, 1, []));   % nq x 1 x N
    uy = pagemtimes(l_e(:, 1:d1, :), reshape(u_hat{2}(:, f), d1, 1, []));
    rhs(ra, 1, :) = len .* pagemtimes(wl, 'transpose', nx .* ux + ny .* uy, 'none');

    % <d_t(u*.n) - t.({L}^T n), d_t mu> = 0, mu the degree-(k+1) mode
    rb = rows_a + e;
    dmu = w1 .* dl_e(:, end, :) ./ len;                 % nq x 1 x N, d_t mu (up to sign)
    dtP1 = tx .* P1ex(nodes_e, :, :) + ty .* P1ey(nodes_e, :, :);   % t . grad Phi
    A(rb, 1:n, :)     = len .* nx .* pagemtimes(dmu, 'transpose', dtP1, 'none');
    A(rb, n+1:2*n, :) = len .* ny .* pagemtimes(dmu, 'transpose', dtP1, 'none');
    % {L} at the face points: own side plus the neighbour across face f
    L_own = evaluate_gradient(Pke(nodes_e, :), L_h, 1:N_ELEMENTS);
    Kn = neighbour_across(mesh, f);
    r_n = pagemldivide(F(:, :, Kn), xy_e(:, nodes_e, :) - reshape(xy_0(:, Kn), 2, 1, []));   % 2 x nq x N
    Pkn = dubiner_basis(k, reshape(r_n(1, :, :), [], 1), reshape(r_n(2, :, :), [], 1));
    Pkn = reshape(Pkn, nq, N_ELEMENTS, d2);
    Pkn = permute(Pkn, [1 3 2]);                        % nq x d2 x N
    L_nb = evaluate_gradient(Pkn, L_h, Kn);
    L_avg = cellfun(@(a, b) (a + b)/2, L_own, L_nb, 'UniformOutput', false);
    % t . L^T n = sum_ij n_i L_ij t_j
    tLn = nx .* (L_avg{1,1} .* tx + L_avg{1,2} .* ty) + ny .* (L_avg{2,1} .* tx + L_avg{2,2} .* ty);
    rhs(rb, 1, :) = len .* pagemtimes(dmu, 'transpose', tLn, 'none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = pagemldivide(A, rhs);
u_star = {reshape(c(1:n, 1, :), n, N_ELEMENTS); reshape(c(n+1:2*n, 1, :), n, N_ELEMENTS)};

if nargout > 1
    div = pagemtimes(P1x, reshape(u_star{1}, n, 1, [])) + pagemtimes(P1y, reshape(u_star{2}, n, 1, []));
    divergence = sqrt(max(vol .* sum(w .* div.^2, 1), [], 'all'));
end
end

function L = evaluate_gradient(P, L_h, K)
% L_h components at the points of P (nq x d2 [x N]) for elements K
L = cell(2, 2);
for i = 1:2
    for j = 1:2
        L{i, j} = pagemtimes(P, reshape(L_h{i, j}(:, K), [], 1, numel(K)));
    end
end
end

function Kn = neighbour_across(mesh, f)
% the element sharing face f with each element of the column (itself on
% boundary faces, where the incidence table repeats the entry)
K = reshape(mesh.incidence(:, 1, f), 2, []);
own = 1:numel(f);
Kn = K(1, :);
Kn(Kn == own) = K(2, Kn == own);
end
