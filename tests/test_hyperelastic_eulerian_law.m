% test_hyperelastic_eulerian_law.m
% The material "neo_hookean_eulerian" of HyperelasticProblem: the Cauchy
% stress sigma_E(G) of the neo-Hookean law as a function of the Eulerian
% displacement gradient G = grad_x u (F = (I - G)^-1), E = 10, nu = 0.3.
%   (a) sigma_E(G) = J^-1 P(F) F' with P the class's Lagrangian law
%       "neo_hookean" (1e-12 relative), for random G with |G| ~ 0.3;
%   (b) sigma_E is symmetric (1e-13) and vanishes on rigid rotations,
%       G = I - Q' (1e-13);
%   (c) its derivative against central differences (1e-6 relative), and
%       D sigma_E(0) = the Hooke tensor (1e-13);
%   (d) the Newton linearization of the class: with Data.iterate = X the
%       local residual newton_residual(..., local_rows = true) at X + eps dX
%       differs from r(X) + eps J dX (J the first local operator of
%       get_matrices) by O(eps^2): the ratio of the defects for eps and
%       eps/2 is 4 (3.5 to 4.5), k = 1, 2.
failures = 0;
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(10, 0.3);
law = @(G) HyperelasticProblem.first_piola(G, lambda, mu, "neo_hookean_eulerian");
mm = @(A, B) [A(:, 1).*B(:, 1) + A(:, 3).*B(:, 2), A(:, 2).*B(:, 1) + A(:, 4).*B(:, 2), ...
              A(:, 1).*B(:, 3) + A(:, 3).*B(:, 4), A(:, 2).*B(:, 3) + A(:, 4).*B(:, 4)];
I = [1 0 0 1];

% (a) identity with the Lagrangian law
rng(7);
G = 0.3*randn(500, 4)/2;
a = I - G;
j = a(:, 1).*a(:, 4) - a(:, 2).*a(:, 3);
G = G(j > 0.2, :);
a = I - G;
j = a(:, 1).*a(:, 4) - a(:, 2).*a(:, 3);
F = [a(:, 4), -a(:, 2), -a(:, 3), a(:, 1)]./j;
P = HyperelasticProblem.first_piola(F - I, lambda, mu, "neo_hookean");
reference = j.*mm(P, F(:, [1 3 2 4]));
sigma = law(G);
err = max(abs(sigma - reference), [], "all")/max(abs(reference), [], "all");
fprintf("(a) %d points, max |G| = %.3f: |sigma_E - J^-1 P F'|/|sigma| = %.2e\n", ...
    size(G, 1), max(abs(G), [], "all"), err);
failures = failures + check(err <= 1e-12, "(a) sigma_E = J^-1 P F'");

% (b) symmetry and rigid rotations
scale = max(abs(sigma), [], "all");
err = max(abs(sigma(:, 2) - sigma(:, 3)))/scale;
fprintf("(b) |sigma_21 - sigma_12|/|sigma| = %.2e\n", err);
failures = failures + check(err <= 1e-13, "(b) symmetry");
theta = linspace(-3, 3, 41).';
Qt = [cos(theta), -sin(theta), sin(theta), cos(theta)];      % Q', column-major
err = max(abs(law(I - Qt)), [], "all")/(2*mu + lambda);
fprintf("(b) rotations: max |sigma_E(I - Q')|/(2 mu + lambda) = %.2e\n", err);
failures = failures + check(err <= 1e-13, "(b) objectivity");

% (c) derivative
[~, dsigma] = law(G);
h = 1e-6;
err = 0;
for b = 1:4
    e = h*((1:4) == b);
    fd = (law(G + e) - law(G - e))/(2*h);
    err = max(err, max(abs(fd - dsigma(:, :, b)), [], "all"));
end
err = err/max(abs(dsigma), [], "all");
fprintf("(c) |dsigma - central differences|/|dsigma| = %.2e\n", err);
failures = failures + check(err <= 1e-6, "(c) derivative");
[~, d0] = law(zeros(1, 4));
hooke = [lambda + 2*mu, 0, 0, lambda; 0, mu, mu, 0; 0, mu, mu, 0; lambda, 0, 0, lambda + 2*mu];
err = max(abs(reshape(d0, 4, 4) - hooke), [], "all")/(lambda + 2*mu);
fprintf("(c) |D sigma_E(0) - C|/(lambda + 2 mu) = %.2e\n", err);
failures = failures + check(err <= 1e-13, "(c) Hooke tensor at G = 0");

% (d) Newton linearization through the class
mesh = StructuredSimplicialGrid(4, 4).add_stabilization(16);
N = size(mesh.elements, 2);
for k = 1:2
    d1 = k + 1;
    d2 = nchoosek(k + 2, 2);
    pde = HyperelasticProblem(k, material = "neo_hookean_eulerian", young_modulus = 10, poisson_ratio = 0.3);
    pde.Order = k;
    [~, ~, rhs_rows] = get_assemble_indices(pde, mesh, k);
    indices = struct("rhs_rows", rhs_rows);
    rng(10 + k);
    skeleton = 1e-3*randn(2*d1, size(mesh.faces, 2));
    u_K = [reshape(skeleton(1:d1, mesh.boundaries), 3*d1, N); ...
           reshape(skeleton(d1 + 1:end, mesh.boundaries), 3*d1, N)];
    x_hat = zeros(max(rhs_rows(:)), 1);
    x_hat(rhs_rows(:)) = u_K(:);
    d2peers = 4*d2 + k + 1;                         % PEERS: 4 Dubiner blocks and k + 1 bubbles
    n_local = d2peers + 3*d2;
    X = 1e-3*randn(n_local, N);
    dX = 1e-3*randn(n_local, N);
    pde.Data.iterate = struct("displacement", X(d2peers + (1:2*d2), :), "displacement_skeleton", skeleton);
    local_operators = get_matrices(pde, mesh, k);
    J = local_operators{1};
    r0 = newton_residual(pde, mesh, k, X, x_hat, indices, local_rows = true);
    defect = zeros(1, 2);
    for i = 1:2
        epsilon = 0.5/2^(i - 1);
        r = newton_residual(pde, mesh, k, X + epsilon*dX, x_hat, indices, local_rows = true);
        linear = r0 + reshape(pagemtimes(J, reshape(epsilon*dX, n_local, 1, N)), n_local, N);
        defect(i) = max(abs(r - linear), [], "all");
    end
    ratio = defect(1)/defect(2);
    fprintf("(d) k = %d: x_hat(rhs_rows) = traces of the skeleton: %d; defects %.3e, %.3e, ratio %.3f\n", ...
        k, isequal(reshape(x_hat(rhs_rows), size(u_K)), u_K), defect, ratio);
    failures = failures + check(isequal(reshape(x_hat(rhs_rows), size(u_K)), u_K), sprintf("(d) k = %d: skeleton layout", k));
    failures = failures + check(ratio >= 3.5 && ratio <= 4.5, sprintf("(d) k = %d: second-order defect", k));
end

% (e) "saint_venant_kirchhoff_eulerian": (a)-(c) for sigma_E(G) = j F S F'
% against the Lagrangian "saint_venant_kirchhoff" (P = F S)
svk = @(G) HyperelasticProblem.first_piola(G, lambda, mu, "saint_venant_kirchhoff_eulerian");
G = 0.3*randn(500, 4)/2;
a = I - G;
j = a(:, 1).*a(:, 4) - a(:, 2).*a(:, 3);
G = G(j > 0.2, :);
a = I - G;
j = a(:, 1).*a(:, 4) - a(:, 2).*a(:, 3);
F = [a(:, 4), -a(:, 2), -a(:, 3), a(:, 1)]./j;
P = HyperelasticProblem.first_piola(F - I, lambda, mu, "saint_venant_kirchhoff");
reference = j.*mm(P, F(:, [1 3 2 4]));
sigma = svk(G);
err = max(abs(sigma - reference), [], "all")/max(abs(reference), [], "all");
fprintf("(e) SVK: |sigma_E - J^-1 P F'|/|sigma| = %.2e\n", err);
failures = failures + check(err <= 1e-12, "(e) SVK sigma_E = J^-1 P F'");
scale = max(abs(sigma), [], "all");
err = max(abs(sigma(:, 2) - sigma(:, 3)))/scale;
err_rot = max(abs(svk(I - Qt)), [], "all")/(2*mu + lambda);
fprintf("(e) SVK: symmetry %.2e, rotations %.2e\n", err, err_rot);
failures = failures + check(err <= 1e-13 && err_rot <= 1e-13, "(e) SVK symmetry and objectivity");
[~, dsigma] = svk(G);
err = 0;
for b = 1:4
    e = h*((1:4) == b);
    fd = (svk(G + e) - svk(G - e))/(2*h);
    err = max(err, max(abs(fd - dsigma(:, :, b)), [], "all"));
end
err = err/max(abs(dsigma), [], "all");
[~, d0] = svk(zeros(1, 4));
err0 = max(abs(reshape(d0, 4, 4) - hooke), [], "all")/(lambda + 2*mu);
fprintf("(e) SVK: |dsigma - central differences|/|dsigma| = %.2e, |D sigma_E(0) - C| = %.2e\n", err, err0);
failures = failures + check(err <= 1e-6 && err0 <= 1e-13, "(e) SVK derivative");

if failures == 0
    fprintf("test_hyperelastic_eulerian_law: all checks passed\n");
else
    fprintf("test_hyperelastic_eulerian_law: %d check(s) failed\n", failures);
end
