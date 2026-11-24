function C = cockburn_shen_constants(mesh, k)
% The local stabilization constants C_KF,o of Cockburn and Shen (Results
% in Applied Mathematics 1 (2019) 100001, Sec. 2.5.5), d = 2:
%   C_KF,o = (k+2)/(k+1) V_F,K + h_F (k+1)(k+2)(k+d)/d |∂K|/|K|,
%   V_F,K  = η_K,F N(k,d) h_F/|F|,   η_K,F = max{ω_V(K), Q(k,d,K)}
% (ω_E enters only in 3D), with
%   Q(k,d,K) the largest eigenvalue of the Lagrange stiffness matrix
%            (∇ψ_i, ∇ψ_j)_K of P_k(K) (equispaced nodes),
%   N(k,d)   of their Table 1 (k = 1, 2, 3),
%   ω_V(K)   = max over the vertices v of F of ω_v,
%   ω_v      = Σ_{K'∈K_v} Q(k,d,K') W_v(K'),
%   W_v(K)   = Σ_{K'∈K_v} Σ_{F'∈F_v∩F(K')} (d_max(K,v) − d_v(K',K) + 1)² / (card F_v)²,
% K_v, F_v the elements and faces containing the vertex v, d_v the
% distance in the graph of K_v whose edges are shared faces, d_max(K,v)
% the largest distance from K (Sec. 2.5.3). Returns C, 3 x N (local face
% e of element K), to be divided by h_F = |F|.
N_table = [6.00, 12.43, 20.77, 30.97];
assert(k >= 1 && k <= numel(N_table), "HyperelasticGradientProblem:degree", ...
    "The Cockburn-Shen constants are tabulated for k = 1, ..., 4.");
N_kd = N_table(k);
N_ELEMENTS = size(mesh.elements, 2);
N_VERTICES = size(mesh.coordinates, 2);

Q = lagrange_stiffness_max(mesh, k);                  % 1 x N
face_length = mesh.area(mesh.boundaries);             % 3 x N
perimeter = sum(face_length, 1);

% local face e of K is opposite its vertex e: its endpoints are the other two
omega_v = zeros(1, N_VERTICES);
element_lists = accumarray(mesh.elements(:), reshape(repmat(1:N_ELEMENTS, 3, 1), [], 1), ...
    [N_VERTICES, 1], @(x) {x});
for v = 1:N_VERTICES
    Kv = element_lists{v};
    if isempty(Kv), continue, end
    n = numel(Kv);
    % faces of each element in K_v containing v (two per element)
    faces = zeros(n, 2);
    for i = 1:n
        e_v = find(mesh.elements(:, Kv(i)) == v);
        faces(i, :) = mesh.boundaries(setdiff(1:3, e_v), Kv(i));
    end
    n_faces = numel(unique(faces));
    % graph: elements adjacent if they share a face containing v
    A = false(n);
    for i = 1:n
        for j = i+1:n
            A(i, j) = any(ismember(faces(i, :), faces(j, :)));
        end
    end
    D = distances(graph(A | A.'));
    W = sum((max(D, [], 2) - D + 1).^2, 2).'*2/n_faces^2;  % two faces of F_v per element
    omega_v(v) = sum(Q(Kv).*W);
end

C = zeros(3, N_ELEMENTS);
for e = 1:3
    ends = mesh.elements(setdiff(1:3, e), :);              % 2 x N endpoints of face e
    omega_V = max(omega_v(ends), [], 1);
    eta = max(omega_V, Q);
    V = eta*N_kd;                                          % h_F/|F| = 1
    C(e, :) = (k + 2)/(k + 1)*V + face_length(e, :)*(k + 1)*(k + 2)*(k + 2)/2 .* perimeter./mesh.volume(:).';
end
end

function Q = lagrange_stiffness_max(mesh, k)
% the largest eigenvalue of (∇ψ_i, ∇ψ_j)_K, ψ the Lagrange basis of P_k(K)
% on the equispaced nodes, per element
[r, s] = meshgrid(0:k);
keep = r(:) + s(:) <= k;
nodes = [r(keep), s(keep)]/k;
V = dubiner_basis(k, nodes(:, 1), nodes(:, 2));          % V_ij = p_j(x_i)
[xq, wq] = get_2d_reference_quadrature(k + 1);
[~, pr, ps] = dubiner_basis(k, xq(1, :).', xq(2, :).');
Lr = pr/V; Ls = ps/V;                                    % reference gradients of ψ
Ginv = get_affine_map(mesh, inverse=true);
N_ELEMENTS = size(mesh.elements, 2);
Q = zeros(1, N_ELEMENTS);
for K = 1:N_ELEMENTS
    gx = Lr*Ginv(1, 1, K) + Ls*Ginv(2, 1, K);
    gy = Lr*Ginv(1, 2, K) + Ls*Ginv(2, 2, K);
    S = mesh.volume(K)*(gx.'*(wq(:).*gx) + gy.'*(wq(:).*gy));
    Q(K) = max(eig((S + S.')/2));
end
end
