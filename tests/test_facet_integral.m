% test_facet_integral.m
% Regression test for the facet and line integral machinery, against the
% boundary_matrix / mass_matrix builders in the elemental page convention
% [d, d, 1, 1, el(, 3)]. Tolerance 1e-12 throughout.
%   Tests 1-5: dS (element-boundary) forms for the vol-vol, hyb-vol,
%              vol-hyb and hyb-hyb trial/test pairs, and a Datum linear form.
%   Tests 6-7: ds (skeleton) forms: P0 facet mass equals |e|; the two
%              element-side copies of each face matrix agree and equal
%              the face mass matrix.
%   Test 8:    LineIntegral.evaluate reproduces I(x) = int_[x, x+d] u dsigma
%              exactly for u = x^2 + y on every face.
clc; clear;

k = 2;
mesh = StructuredSimplicialGrid(2, 2);

Pk = DubinerBasis(k);
Mk = LegendreBasis(k);

% Trusted reference matrices; PP: [n n 1 1 el 3], LP: [n1 n 1 1 el 3]
face_pages = @(A) permute(A, [1 2 5 6 4 3]);
pages = @(A) permute(A, [1 2 4 5 3]);
PP = face_pages(boundary_matrix(mesh, Pk, Pk));
LP = face_pages(boundary_matrix(mesh, Pk, Mk));
LL = pages(mass_matrix(mesh, Mk, Mk));

dS = Measure(mesh, 'dS');

n_fail = 0;

%% Test 1: volume-volume boundary mass  sum_e int_e phi_i phi_j
u = TrialFunction(Pk);
w = TestFunction(Pk);
try
    A = double(u * w * dS);
    A_ref = sum(PP, 6);
    err = max(abs(A(:) - A_ref(:)));
    fprintf("Test 1 (u*w*dS, vol-vol): size = [%s], err = %.3e\n", ...
        num2str(size(A)), err);
    if err > 1e-12; n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 1 (u*w*dS) ERROR: %s\n", ME.message);
end

%% Test 2: hybrid-volume boundary mass  int_e lambda_i phi_j
U = TrialFunction(Mk);
try
    B = double(U * w * dS);
    % reference: LP has rows = Legendre, cols = Dubiner; trial U indexes cols
    % of the local matrix (I(i,j) = int lambda_j phi_i) => transpose of LP
    B_ref = permute(pagetranspose(LP), [1 2 3 4 5 7 6]);
    fprintf("Test 2 (U*w*dS, hyb-vol): size B = [%s], size ref = [%s]\n", ...
        num2str(size(B)), num2str(size(B_ref)));
    err = max(abs(B(:) - B_ref(:)));
    fprintf("        err = %.3e\n", err);
    if err > 1e-12; n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 2 (U*w*dS) ERROR: %s\n", ME.message);
end

%% Test 3: volume-hybrid boundary mass  int_e phi_j lambda_i
V = TestFunction(Mk);
try
    C = double(u * V * dS);
    C_ref = LP;
    fprintf("Test 3 (u*V*dS, vol-hyb): size C = [%s], size ref = [%s]\n", ...
        num2str(size(C)), num2str(size(C_ref)));
    err = max(abs(C(:) - C_ref(:)));
    fprintf("        err = %.3e\n", err);
    if err > 1e-12; n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 3 (u*V*dS) ERROR: %s\n", ME.message);
end

%% Test 4: hybrid-hybrid boundary mass  int_e lambda_i lambda_j
try
    D = double(U * V * dS);
    D_ref = reshape(LL(:, :, :, :, mesh.boundaries.'), ...
                    size(LL, 1), size(LL, 2), 1, 1, [], 3);
    fprintf("Test 4 (U*V*dS, hyb-hyb): size D = [%s], size ref = [%s]\n", ...
        num2str(size(D)), num2str(size(D_ref)));
    % note ref block-diagonal structure expected in dims 6-7
    err = NaN;
    if isequal(size(D, 1:5), size(D_ref, 1:5)) && size(D, 6) == 3 && size(D, 7) == 3
        D_diag = zeros(size(D_ref));
        for e = 1:3; D_diag(:, :, :, :, :, e) = D(:, :, :, :, :, e, e); end
        err = max(abs(D_diag(:) - D_ref(:)));
    end
    fprintf("        err (diag blocks) = %.3e\n", err);
    if ~(err < 1e-12); n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 4 (U*V*dS) ERROR: %s\n", ME.message);
end

%% Test 5: linear form  int_dK g w dS with g = 1 (exact: int of basis over dK)
try
    g = Datum(@(x, y) 1 + 0*x);
    L = double(g * w * dS);
    fprintf("Test 5 (g*w*dS, datum-vol): size L = [%s]\n", num2str(size(L)));
    % reference: sum_e |e| * int_0^1 phi(edge(s)) ds per element
    [s_local, wq] = get_1d_reference_quadrature(2*k);
    x_edges = [[0*s_local; 1 - s_local], [1 - s_local; 0*s_local], [1 - s_local; s_local]];
    p = Pk.evaluate(x_edges);                          % (3*nq) x n
    p = reshape(double(p), length(s_local), 3, []);    % nq x 3 x n
    a_e = permute(mesh.area(mesh.boundaries), [3 1 4 2]);  % 1 x 3 x 1 x el
    L_ref = sum(sum(wq(:) .* p .* a_e, 1), 2);         % 1 x 1 x n x el
    L_ref = squeeze(L_ref);                            % n x el
    Lsq = squeeze(L);
    err = max(abs(Lsq(:) - L_ref(:)));
    fprintf("        err = %.3e\n", err);
    if ~(err < 1e-12); n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 5 (g*w*dS) ERROR: %s\n", ME.message);
end

%% Test 6: skeleton measure 'ds' — P0 facet mass: int_e 1*1 ds = |e| per face
try
    ds = Measure(mesh, 'ds');
    U0 = TrialFunction(LegendreBasis(0));
    V0 = TestFunction(LegendreBasis(0));
    E = double(U0 * V0 * ds);
    fprintf("Test 6 (U0*V0*ds, skeleton): size E = [%s]\n", num2str(size(E)));
    Esq = squeeze(E);          % n_faces x 2 (one copy per side of the face)
    err = max(abs(Esq - mesh.area(:)), [], 'all');
    fprintf("        err vs |e| = %.3e\n", err);
    if ~(err < 1e-12); n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 6 (U0*V0*ds) ERROR: %s\n", ME.message);
end

%% Test 7: skeleton measure 'ds' with degree-k trace functions
try
    Ek = double(TrialFunction(Mk) * TestFunction(Mk) * ds);
    fprintf("Test 7 (Uk*Vk*ds): size = [%s]\n", num2str(size(Ek)));
    % both sides of an interior face carry the same trace space; the two
    % copies (dim 6) must agree, and each must equal |e| * 1D mass matrix
    LL1 = pages(mass_matrix(mesh, Mk, Mk));
    err_sides = max(abs(Ek(:,:,:,:,:,1) - Ek(:,:,:,:,:,2)), [], 'all');
    err_mass = max(abs(Ek(:,:,:,:,:,1) - LL1), [], 'all');
    fprintf("        max |side1 - side2| = %.3e, max |side1 - |e| M1d| = %.3e\n", ...
        err_sides, err_mass);
    if ~(err_sides < 1e-12 && err_mass < 1e-12); n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 7 (Uk*Vk*ds) ERROR: %s\n", ME.message);
end

%% Test 8: LineIntegral (TPM transfer path) against exact integral
% u(x,y) = x^2 + y (in the FE space for k >= 2), transfer path from each
% face point x to x_bar = x + d with constant offset d. Exact:
%   I(x) = |d| * int_0^1 u(x + s*d) ds
%        = |d| * (x1^2 + x1*d1 + d1^2/3 + x2 + d2/2)
try
    uh = TrialFunction(Pk);
    uh.mesh = mesh;
    % project u(x,y) = x^2 + y onto the (elementwise) Dubiner space:
    % dofs from L2 projection using the trusted mass builder
    M = pages(mass_matrix(mesh, Pk, Pk));
    fh = Datum(@(x, y) x.^2 + y);
    w2 = TestFunction(Pk);
    dx = Measure(mesh, 'dx');
    b = double(fh * w2 * dx);
    b = reshape(b, size(b, 1), 1, []);
    dofs = pagemldivide(M, b);
    field = FEField(Pk, squeeze(dofs));
    field.mesh = mesh;

    d = [0.3; -0.2];
    li = LineIntegral(field, @(X) X + d);
    x_loc = [0.25 0.75];
    out = li.evaluate(x_loc, mesh);
    fprintf("Test 8 (LineIntegral.evaluate): size = [%s]\n", num2str(size(out)));

    x_f = get_physical_points(mesh, x_loc);   % 2 x 2 x n_faces
    x1 = squeeze(x_f(1, :, :)); x2 = squeeze(x_f(2, :, :));  % 2 x n_faces
    I_exact = norm(d) * (x1.^2 + x1*d(1) + d(1)^2/3 + x2 + d(2)/2);
    out_sq = reshape(out, length(x_loc), []);
    err = max(abs(out_sq - I_exact), [], 'all');
    fprintf("        err vs exact = %.3e\n", err);
    if ~(err < 1e-12); n_fail = n_fail + 1; end
catch ME
    n_fail = n_fail + 1;
    fprintf("Test 8 (LineIntegral) ERROR: %s\n", ME.message);
    fprintf("%s\n", getReport(ME, 'basic', 'hyperlinks', 'off'));
end

fprintf("\n%d check(s) failed.\n", n_fail);
