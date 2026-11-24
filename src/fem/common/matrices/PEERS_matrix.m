function [M_PEERS, M_tr_PEERS, M_A_PEERS] = PEERS_matrix(mesh, k)
%PEERS_MATRIX  Element matrices of the PEERS basis.
%
%   [M_PEERS, M_tr_PEERS, M_A_PEERS] = PEERS_matrix(mesh, k)
%
%   With P the PEERS basis of order k (dimension n_P) and p the Dubiner
%   basis of order k (dimension d2), per element:
%     M_PEERS     int P : P            n_P-by-n_P-by-n_el
%     M_tr_PEERS  int tr(P) tr(P)      n_P-by-n_P-by-n_el
%     M_A_PEERS   int P : (p [0 1; -1 0])   n_P-by-d2-by-n_el
arguments
    mesh SimplicialGrid
    k double
end
% The products of mass_matrix's PEERS branch (rows stacked by stack_rows,
% one page product per element) on a single evaluation of the PEERS basis;
% the antisymmetric operand is element-independent, one page broadcast
% over the elements by pagemtimes.
[xhat, w] = get_2d_reference_quadrature(k + 2);
vol = reshape(mesh.volume, 1, 1, []);

P = PEERS_basis(k, xhat(1, :), xhat(2, :), get_affine_map(mesh, inverse=true));
A_P = tensor_basis(k, xhat(1, :), xhat(2, :), ...
                   PolynomialBasis=@dubiner_basis, ...
                   LinearMaps={[0, 1; -1, 0]});

p = stack_rows(P);
p_tr = stack_rows(P(:, :, 1, 1, :) + P(:, :, 2, 2, :));
p_A = stack_rows(A_P);

w_p = repmat(w(:), size(p, 1)/numel(w), 1) .* p;
M_PEERS = vol .* pagemtimes(w_p, 'transpose', p, 'none');
M_tr_PEERS = vol .* pagemtimes(w(:) .* p_tr, 'transpose', p_tr, 'none');
M_A_PEERS = vol .* pagemtimes(w_p, 'transpose', p_A, 'none');
end

function P = stack_rows(P)
% (n_q x dim x c1 x c2 x n_el) -> (n_q*c1*c2 x dim x n_el), quadrature
% point fastest, so the same weights repeat over every component block.
[n_q, dim, c1, c2, n_el] = size(P, 1:5);
P = reshape(permute(P, [1 3 4 2 5]), n_q*c1*c2, dim, n_el);
end
