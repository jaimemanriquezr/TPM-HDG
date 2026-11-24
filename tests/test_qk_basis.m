%% QkBasis: tensor-product Legendre basis on the reference square
% Checks, for k = 0:4: Dimension; agreement with legendre_eval in each
% factor; orthonormality on [0, 1]^2; first and second derivatives
% against centred differences; and the Kronecker structure of the
% reference stiffness matrix, K2 = kron(I, K1) + kron(K1, I).
clc; clear;
n_fail = 0;
tol = 1e-12;
tol_fd = 1e-6;

for k = 0:4
    basis = QkBasis(k);
    m = k + 1;
    check("Dimension", abs(basis.Dimension - m^2), 0);

    % factor agreement: phi_(i,1) = L_i(r) * L_1(s) and L_1 = 1
    [x, w] = get_2d_reference_quadrature(k + 1, "Quad");
    [P, dPx, dPy] = basis.evaluate(x);
    L = legendre_eval(1:m, x(1, :).');
    check("factor", max(abs(P(:, 1:m) - L), [], 'all'), tol);

    % orthonormality: reference mass matrix is the identity
    M = P.' * (w .* P);
    check("mass", max(abs(M - eye(m^2)), [], 'all'), tol);

    % Kronecker structure of the stiffness matrix
    [s1, w1] = get_1d_reference_quadrature(k + 1, "Quad");
    [L1, dL1] = legendre_eval(1:m, s1.');
    K1 = dL1.' * (w1 .* dL1);
    K2 = dPx.' * (w .* dPx) + dPy.' * (w .* dPy);
    K_kron = kron(eye(m), K1) + kron(K1, eye(m));
    check("stiffness kron", max(abs(K2 - K_kron), [], 'all'), tol*max(1, norm(K1, 'inf')));

    % derivatives against centred differences at interior points
    y = [0.13 0.41 0.77; 0.29 0.58 0.91];
    h = 1e-5; ex = [h; 0]; ey = [0; h];
    D = cell(1, 6); [D{:}] = QkBasis.evaluate_basis(k, y);
    Dpx = cell(1, 3); [Dpx{:}] = QkBasis.evaluate_basis(k, y + ex);
    Dmx = cell(1, 3); [Dmx{:}] = QkBasis.evaluate_basis(k, y - ex);
    Dpy = cell(1, 3); [Dpy{:}] = QkBasis.evaluate_basis(k, y + ey);
    Dmy = cell(1, 3); [Dmy{:}] = QkBasis.evaluate_basis(k, y - ey);
    cd_x = @(i) (Dpx{i} - Dmx{i}) / (2*h);   % centred d/dx of output i
    cd_y = @(i) (Dpy{i} - Dmy{i}) / (2*h);
    scale = max(1, max(abs(D{4}), [], 'all'));
    rel = @(a, b) max(abs(a - b), [], 'all') / scale;
    check("dPx", rel(D{2}, cd_x(1)), tol_fd);
    check("dPy", rel(D{3}, cd_y(1)), tol_fd);
    check("dPxx", rel(D{4}, cd_x(2)), tol_fd);
    check("dPyy", rel(D{5}, cd_y(3)), tol_fd);
    check("dPxy", rel(D{6}, cd_y(2)), tol_fd);
end
fprintf("%d check(s) failed\n", n_fail);

function check(name, err, tol)
    k = evalin('caller', 'k');
    ok = err <= tol;
    fprintf("k = %d %-15s err = %.2e  %s\n", k, name, err, string(ok));
    if ~ok
        assignin('caller', 'n_fail', evalin('caller', 'n_fail') + 1);
    end
end
