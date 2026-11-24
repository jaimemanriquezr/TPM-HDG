function [Cb, BbPL, BbLL, max_bn] = convective_matrix(mesh, k, beta, options)
arguments
    mesh SimplicialGrid
    k double
    beta (2,1) cell   % {bx; by}: function handles of (x,y), or d2 x Nelts dof arrays
    options.basis_2d function_handle = @dubiner_basis
    options.basis_1d function_handle = @legendre_basis
    % faces whose convective flux is dropped from the global flux equation
    % (do-nothing outflow: the traction condition is sigma n - S(u - uhat) =
    % g there, not sigma n - rho (beta.n) uhat = g). Their BbLL block is zero;
    % BbPL (the local momentum face term) is kept, so the local equation
    % stays the conservative-form Navier-Stokes momentum balance.
    options.skip_faces (:,1) double = []
end
%-------------------------------------------------------------------------
%Description: Element matrices of the Oseen convective term for a given
%             convective velocity beta (Solano & Vargas 2022, eq. 3.4):
%               Cb_ab   = int_K  P_b (beta . grad P_a)        (volume)
%               BbPL_ca = int_e  (beta . n) L_c P_a           (face/volume)
%               BbLL_cd = int_e  (beta . n) L_c L_d           (face/face)
%             so that, per velocity component,
%               -(u x beta, grad v)_K       ->  -Cb
%                <(uhat x beta) n, v>_dK    ->   BbPL'
%                <(uhat x beta) n, mu>_dK   ->   BbLL
%             beta is either a pair of vectorised handles or a pair of
%             Dubiner dof arrays of degree k (the previous velocity of a
%             Picard or time step). max_bn = max |beta . n| over the face
%             quadrature points, for the stabilization rule
%             tau = max_bn/(2 nu) + 1.
%Output:    - Cb:   (d2 x d2 x Nelts)
%           - BbPL: (3*d1 x d2 x Nelts)
%           - BbLL: (3*d1 x 3*d1 x Nelts), block diagonal per face
%           - max_bn: scalar
%-------------------------------------------------------------------------
N_ELEMENTS = size(mesh.elements, 2);
is_discrete = ~isa(beta{1}, 'function_handle');
kb = k;                       % degree of a discrete beta: P_k, or P_{k+1} (u*)
if is_discrete
    kb = find(arrayfun(@(d) nchoosek(d + 2, 2), 0:k+2) == size(beta{1}, 1)) - 1;
    assert(~isempty(kb), "convective_matrix:betaDegree", "beta dofs match no P_d");
end
kq = max(ceil((kb + 2*k)/2), 1);   % beta in P_kb times two P_k factors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VOLUME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xy, rs, w, B] = quadpoints_elements_2d(mesh, kq, output="cell");
[p, pr, ps] = options.basis_2d(k, rs{:});
d2 = size(p, 2);
[bx, by] = evaluate_beta(beta, xy, options.basis_2d(kb, rs{:}), is_discrete);

% physical gradients: grad P = B^{-T} grad_hat P, det B = 2|K|
det_B = B{1,1}(:).' .* B{2,2}(:).' - B{1,2}(:).' .* B{2,1}(:).';
px = (B{2,2}(:).' .* pr(:) - B{2,1}(:).' .* ps(:));   % (nq*d2) x N, scaled by det
py = (-B{1,2}(:).' .* pr(:) + B{1,1}(:).' .* ps(:));
px = reshape(px, [], d2, N_ELEMENTS) ./ reshape(det_B, 1, 1, []);
py = reshape(py, [], d2, N_ELEMENTS) ./ reshape(det_B, 1, 1, []);

% beta . grad P_a at each quadrature point, then against w P_b
bgp = reshape(bx, [], 1, N_ELEMENTS) .* px + reshape(by, [], 1, N_ELEMENTS) .* py;
Cb = pagemtimes(bgp, 'transpose', w .* p, 'none');          % (a, b)
Cb = Cb .* reshape(mesh.volume, 1, 1, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s, w1] = get_1d_reference_quadrature(kq);
N_NODES = length(w1);
quadrature_2d = [[0*s(:), 1 - s(:)];
                 [1 - s(:), 0*s(:)];
                 [1 - s(:), s(:)]];
pe = options.basis_2d(k, quadrature_2d(:,1), quadrature_2d(:,2));  % 3*nq x d2
l = options.basis_1d(k, [s(:); 1 - s(:)]);                          % 2*nq x d1
d1 = size(l, 2);

[F, xy_0] = get_affine_map(mesh, output="array");
xy_e = pagemtimes(F, quadrature_2d.') + reshape(xy_0, 2, 1, []);   % 2 x 3nq x N
xy_e = {reshape(xy_e(1, :, :), [], N_ELEMENTS); reshape(xy_e(2, :, :), [], N_ELEMENTS)};
pb_e = options.basis_2d(kb, quadrature_2d(:,1), quadrature_2d(:,2));
[bx_e, by_e] = evaluate_beta(beta, xy_e, pb_e, is_discrete);        % 3nq x N

% mesh.normals carry the face length, so no area factor (as in nxLP)
area = mesh.area(mesh.boundaries);
oriented_normals = mesh.normals .* kron(mesh.orientation, [1;1]);
BbPL = zeros(3*d1, d2, N_ELEMENTS);
BbLL = zeros(3*d1, 3*d1, N_ELEMENTS);
max_bn = 0;
for e = 1:3
    nodes_e = (e-1)*N_NODES + (1:N_NODES);
    block_e = (e-1)*d1 + (1:d1);
    bn = oriented_normals(2*(e-1) + 1, :) .* bx_e(nodes_e, :) ...
       + oriented_normals(2*(e-1) + 2, :) .* by_e(nodes_e, :);    % nq x N
    max_bn = max(max_bn, max(abs(bn ./ area(e, :)), [], 'all'));
    wbn = reshape(w1 .* bn, N_NODES, 1, N_ELEMENTS);
    keep = reshape(~ismember(mesh.boundaries(e, :), options.skip_faces), 1, 1, []);
    for mu = 1:2
        is_valid_edge = reshape(mesh.permutation(e, :) == mu, 1, 1, []);
        l_mu = l((mu - 1)*N_NODES + (1:N_NODES), :);
        wl = wbn .* is_valid_edge .* l_mu;                          % nq x d1 x N
        BbPL(block_e, :, :) = BbPL(block_e, :, :) + pagemtimes(wl, 'transpose', pe(nodes_e, :), 'none');
        BbLL(block_e, block_e, :) = BbLL(block_e, block_e, :) + pagemtimes(wl .* keep, 'transpose', l_mu, 'none');
    end
end
end

function [bx, by] = evaluate_beta(beta, xy, p, is_discrete)
% beta at quadrature points: (n_points x Nelts) for each component
if is_discrete
    bx = p * beta{1};
    by = p * beta{2};
else
    bx = beta{1}(xy{:});
    by = beta{2}(xy{:});
    bx = bx + 0*xy{1};  by = by + 0*xy{1};   % constants broadcast
end
end
