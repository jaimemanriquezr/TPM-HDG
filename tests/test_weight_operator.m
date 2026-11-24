% test_weight_operator.m
% WeightOperator forms: constant scalar/matrix, spatial (Datum) scalar and
% matrix-valued weights, and the scalar-operand x matrix-weight promotion.
%
% Check 1: constant matrix weight == the same matrix as constant handles
%          through the spatial path (exact, quadrature-independent).
% Check 2: variable scalar weight vs manual quadrature reference (same rule).
% Check 3: matrix-valued variable weight vs blockwise scalar references.
% Check 4: promotion: matrix weight on a scalar basis gives A(i,j)*phi.
clc; clear;

k = 2;
mesh = StructuredSimplicialGrid(4, 4);
Pk = DubinerBasis(k);
VPk = Pk.to_vector();
n_el = size(mesh.elements, 2);
n_failed = 0;

%% Check 1: constant matrix, double vs spatial-handle path
K_mat = [2 1; 0 3];
M_const = mass_matrix(mesh, VPk, trial_operator = WeightOperator(K_mat));
K_cell = {@(x,y) 2+0*x, @(x,y) 1+0*x; @(x,y) 0*x, @(x,y) 3+0*x};
M_spatial = mass_matrix(mesh, VPk, trial_operator = WeightOperator(K_cell));
err1 = max(abs(M_const - M_spatial), [], 'all');
fprintf("Check 1 (constant matrix, double vs spatial path): err = %.3e\n", err1);
n_failed = n_failed + (err1 > 1e-13);

%% Check 2: variable scalar weight vs manual quadrature reference
w_fun = @(x, y) 1 + x.^2 + 2*y;
M_w = mass_matrix(mesh, Pk, trial_operator = WeightOperator(w_fun));

[x_local, w_q] = get_quadrature([Pk, Pk]);
phi = double(Pk.evaluate(x_local));                 % (n_q x d)
xq = get_2d_physical(mesh, x_local);                % (2 x n_q x n_el)
M_ref = zeros(size(M_w));
for K = 1:n_el
    w_vals = w_fun(xq(1, :, K), xq(2, :, K)).';     % (n_q x 1)
    M_ref(:, :, K) = mesh.volume(K) * (phi.' * (w_q(:) .* w_vals .* phi));
end
err2 = max(abs(M_w - M_ref), [], 'all');
fprintf("Check 2 (variable scalar weight vs quadrature ref): err = %.3e\n", err2);
n_failed = n_failed + (err2 > 1e-13);

%% Check 3: matrix-valued variable weight vs blockwise scalar references
kappa_cell = {@(x,y) 1+x.^2, @(x,y) x.*y; @(x,y) 0*x, @(x,y) 1+y.^2};
M_K = mass_matrix(mesh, VPk, trial_operator = WeightOperator(kappa_cell));
d = Pk.Dimension;
M_blocks = zeros(size(M_K));
for a = 1:2
    for b = 1:2
        M_ab = mass_matrix(mesh, Pk, ...
            trial_operator = WeightOperator(kappa_cell{a, b}));
        M_blocks((a-1)*d + (1:d), (b-1)*d + (1:d), :) = M_ab;
    end
end
err3 = max(abs(M_K - M_blocks), [], 'all');
fprintf("Check 3 (matrix variable weight vs scalar blocks):  err = %.3e\n", err3);
n_failed = n_failed + (err3 > 1e-13);

%% Check 4: scalar-operand x matrix-weight promotion
A = [0 1; -1 0];
W = WeightOperator(A);
P = W.evaluate(Pk, x_local, mesh);                  % (n_q x d x 2 x 2)
err4 = 0;
for i = 1:2
    for j = 1:2
        err4 = max(err4, max(abs(P(:, :, i, j) - A(i, j)*phi), [], 'all'));
    end
end
fprintf("Check 4 (scalar x matrix promotion):               err = %.3e\n", err4);
n_failed = n_failed + (err4 > 1e-13);

fprintf("\n%d check(s) failed.\n", n_failed);
