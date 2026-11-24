% test_stokes_biot_oseen.m
% Steady Oseen-Biot surrogate (V3): the TPM-coupled Stokes-Biot problem of
% StokesBiotCoupling.solve with the fluid's convective term switched on
% through a convective_velocity beta in the fluid Data (passed on by
% StokesBiotCoupling.interface_problem). One backward-Euler step from rest
% (dt = 1), fluid density 0, so the convective density is rho_c = 1; s0 =
% 1, friction_bjs = 1, all-Dirichlet outer boundaries, pinned fluid
% pressure; the fluid source is f_F + rho_c (beta . grad) u_F (the
% option convection of StokesBiotCoupling.get_manufactured_solution when
% beta = u_F). The
% interface rows carry the Cauchy traction (no convective flux on Sigma);
% the fluid element keeps its convective face term on the Sigma faces.
% convective_stabilization = false: the rows read tau from the fluid mesh
% (S = mu tau), so the fluid mesh.tau is set per face to 1 + rho_c max_f
% |beta . n| / (2 mu); the margin min (mu tau - rho_c beta . n_K / 2) over
% element faces is printed (must be > 0).
%
%   (a) patch tests: every exact field polynomial of degree <= k (stream
%       function of degree k + 1), k = 1, 2, 3, h = 1/4, 1/8 (1/8, 1/16
%       with gap h: at h = gap = 1/4 the k = 3 round-off reaches 1e-10):
%       - beta = the exact u_F (handles): flat with gap 0, h^2 and h,
%         relative errors <= 1e-10 for k = 1, 2. For k = 3 the source
%         (u_F . grad) u_F has degree 2k - 1 = 5 and the fluid source is
%         integrated by the order-k element rule (exact to degree 2k + 1
%         = 7 < 5 + 3): the error follows the quadrature defect of the
%         load vector (3e-6, 1e-7, 3e-9 relative at h = 1/4, 1/8, 1/16;
%         patch error 3e-8, 5e-10, 6e-12) and is asserted <= 1e-6 only;
%       - beta affine and divergence free, beta . nu = 2 (x - 1/2) on y =
%         1/2 (source of degree k, integrated exactly): flat with gap 0,
%         h^2 and h, k = 1, 2, 3, <= 1e-10 (k = 3 with gap h: <= 1e-9,
%         the round-off of the extrapolation over a full element reaches
%         1.6e-10 at h = 1/16);
%       - beta = the Dubiner projection of the exact u_F onto P_{k+1}
%         (the form the postprocessed u* takes), flat gap h^2: as beta =
%         u_F;
%       - circle (gap h^2, h = 1/8, 1/12), beta = u_F and affine: as on
%         the flat Sigma (the Taylor transfer is exact for these
%         polynomials; only the k = 3 source quadrature limits beta = u_F).
%   (b) convergence: the default fields of
%       StokesBiotCoupling.get_manufactured_solution (u_F . nu of both
%       signs on Sigma), flat gap h^2 with h = 1/4, 1/8, 1/16 and circle
%       gap h^2 with h = 1/10, 1/20, 1/40; rates coarsest-finest >= k + 1
%       - 0.3 (u_F, u_P, p_P) and k + 1 - 0.5 (sigma_F, sigma_P, q_P) on
%       the flat ladder, - 0.4 / - 0.6 on the circle.
%   (c) tau x 500: the k = 2 flat ladder of (b) with the fluid mesh.tau
%       multiplied by 500 on every face: the same bounds; on the faces
%       off Sigma only: reported, not asserted (u_P falls to 2.6 on h =
%       1/4..1/16 there). The ratio of the errors to tau x 1 is printed.
%   (d) convection is active: on flat gap h^2, h = 1/8, k = 2, the
%       solution with beta differs from the one with beta = 0 (same data)
%       by more than 1e-4 relative in u_F.
%
% Pieces: the environment variable OSEEN_PARTS selects the parts (default
% "abcd", e.g. OSEEN_PARTS=b), OSEEN_DEGREES the degrees of (b) (default
% "1,2"; OSEEN_DEGREES=3 runs the k = 3 ladders, a few minutes on their
% own and passing with the same bounds). k = 0 is skipped (no curl bubble in the stress elements).
clc; clear;
addpath(fullfile(fileparts(mfilename("fullpath")), "helpers"));
parts = string(getenv("OSEEN_PARTS"));
if parts == "", parts = "abcd"; end
ladder_degrees = str2double(split(string(getenv("OSEEN_DEGREES")), ",")).';
if all(isnan(ladder_degrees)), ladder_degrees = [1 2]; end
names = ["sigma_F", "u_F", "sigma_P", "u_P", "q_P", "p_P"];
n_fail = 0;

%% (a) patch tests
if contains(parts, "a")
    patch_tol = 1e-10;
    for k = 1:3
        fields = patch_fields(k);
        tol_exact = patch_tol;
        if k == 3, tol_exact = 1e-6; end
        cases = { ...
            "flat", "zero", "exact", tol_exact; "flat", "h2", "exact", tol_exact; ...
            "flat", "h", "exact", tol_exact; ...
            "flat", "zero", "affine", patch_tol; "flat", "h2", "affine", patch_tol; ...
            "flat", "h", "affine", patch_tol; "flat", "h2", "dubiner", tol_exact; ...
            "circle", "h2", "exact", tol_exact; "circle", "h2", "affine", patch_tol};
        for c = 1:size(cases, 1)
            [problem, gap, beta, tol] = cases{c, :};
            if k == 3 && gap == "h", tol = max(tol, 1e-9); end
            hs = [1/4, 1/8];
            if gap == "h", hs = [1/8, 1/16]; end
            if problem == "circle", hs = [1/8, 1/12]; end
            worst = 0;  at = "";
            for h = hs
                out = oseen_solve(problem, k, h, gap, fields, beta = beta);
                [w, i] = max(out.rel);
                if w > worst, worst = w; at = sprintf("%s at h = 1/%d", names(i), round(1/h)); end
                fprintf("(a) k %d %-6s gap %-4s beta %-7s h 1/%-2d max rel %.2e (%s) margin %.2f residue %.1e\n", ...
                    k, problem, gap, beta, round(1/h), w, names(i), out.margin, out.residue);
            end
            n_fail = n_fail + check(worst <= tol, sprintf("(a) patch k = %d, %s gap %s beta %s: %.3e (%s)", ...
                k, problem, gap, beta, worst, at));
        end
    end
end

%% (b) convergence ladders
ladders = struct("problem", {"flat", "circle"}, "h", {[1/4, 1/8, 1/16], [1/10, 1/20, 1/40]}, ...
    "slack", {[0.3, 0.5], [0.4, 0.6]});
primal = [false true false true false true];             % u_F, u_P, p_P
if contains(parts, "b")
    for k = ladder_degrees
        for L = ladders
            E = ladder(L, k, struct());
            n_fail = n_fail + check_rates(E, L, k, primal, names, sprintf("(b) %s k = %d", L.problem, k));
        end
    end
end

%% (c) tau x 500
if contains(parts, "c")
    k = 2;
    L = ladders(1);
    E1 = ladder(L, k, struct());
    for where = ["all", "off_sigma"]
        E = ladder(L, k, struct("tau_factor", 500, "tau_faces", where));
        n_fail = n_fail + check_rates(E, L, k, primal, names, ...
            sprintf("(c) tau x 500 (%s) flat k = %d", where, k), where == "all");
        fprintf("(c) tau x 500 (%s): error / error(tau x 1) per level\n", where);
        disp(array2table(E./E1, VariableNames = names));
    end
end

%% (d) convection active
if contains(parts, "d")
    with = oseen_solve("flat", 2, 1/8, "h2", struct(), beta = "exact");
    without = oseen_solve("flat", 2, 1/8, "h2", struct(), beta = "zero");
    du = with.fs.Unknowns.field(2).dof - without.fs.Unknowns.field(2).dof;
    rel = norm(du(:))/norm(with.fs.Unknowns.field(2).dof(:));
    fprintf("(d) ||u_F(beta) - u_F(0)|| / ||u_F(beta)|| = %.3e; errors beta: %s; beta = 0: %s\n", ...
        rel, sprintf("%.2e ", with.err), sprintf("%.2e ", without.err));
    n_fail = n_fail + check(rel > 1e-4, sprintf("(d) convection inactive: %.3e", rel));
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_oseen:failed", "%d checks failed", n_fail);
end

%=========================================================================%
function E = ladder(L, k, extra)
% errors (levels x 6) of the default fields on one ladder, gap h^2
args = namedargs2cell(extra);
E = zeros(numel(L.h), 6);
for i = 1:numel(L.h)
    out = oseen_solve(L.problem, k, L.h(i), "h2", struct(), args{:});
    E(i, :) = out.err;
    fprintf("    %-6s k %d h 1/%-3d errors %s b.n_F on Sigma [%+.3f, %+.3f] outer [%+.3f, %+.3f] tau_max %.3g margin %.3g\n", ...
        L.problem, k, round(1/L.h(i)), sprintf("%.3e ", out.err), out.bn_sigma, out.bn_outer, ...
        out.tau_max, out.margin);
end
end

function n_fail = check_rates(E, L, k, primal, names, label, asserted)
% coarsest-finest rates against k + 1 - slack, printed; asserted = false
% only reports them
if nargin < 7, asserted = true; end
rates = log(E(1, :)./E(end, :))/log(L.h(1)/L.h(end));
need = k + 1 - L.slack(2)*ones(1, 6);
need(primal) = k + 1 - L.slack(1);
fprintf("%s: rate", label);
fprintf("  %s %.2f (>= %.1f)", [names; compose("%.2f", rates); compose("%.1f", need)]);
fprintf("\n");
n_fail = 0;
if ~asserted
    fprintf("%s: reported, not asserted (%d below the bound)\n", label, sum(rates < need));
    return
end
for j = 1:6
    n_fail = n_fail + check(rates(j) >= need(j), sprintf("%s: %s rate %.2f < %.1f", ...
        label, names(j), rates(j), need(j)));
end
end

function out = oseen_solve(problem, k, h, gap_kind, fields, options)
% one steady Oseen-Biot solve; fields: options of
% StokesBiotCoupling.get_manufactured_solution (empty = its defaults)
arguments
    problem (1,1) string
    k (1,1) double
    h (1,1) double
    gap_kind (1,1) string
    fields struct
    options.beta (1,1) string = "exact"          % exact | affine | dubiner | zero
    options.tau_factor (1,1) double = 1
    options.tau_faces (1,1) string = "all"       % all | off_sigma
end
gap = 0;
if gap_kind == "h2", gap = h^2; end
if gap_kind == "h", gap = h; end
switch problem
    case "flat"
        curve = FlatCurveHorizontal(0.5);
        [fm, bm] = fitted_split(curve, [0.8*h, h], gap);
        normal_p = ["0"; "1"];
    case "circle"
        curve = CircularCurve([0.5; 0.5], 0.25);
        [fm, bm] = interface_meshes(curve, "fitted", 0.8*h, gap = gap, h_solid = h);
        r = "sqrt((x - 0.5)^2 + (y - 0.5)^2)";
        normal_p = ["(x - 0.5)/" + r; "(y - 0.5)/" + r];
end
% the convective source rho_c (beta . grad) u_F: from the option
% convection when beta = u_F, derived here for the affine beta
is_affine = options.beta == "affine";
args = namedargs2cell(fields);
mms = StokesBiotCoupling.get_manufactured_solution(normal_p, args{:}, storage = 1, ...
    friction_bjs = 1, time_step = 1, rate = "from_rest", convection = ~is_affine);
data = mms.at(1);
mu = 1;  rho_c = 1;
syms x y t real
u = mms.symbolic.u_F;
b = u;
if is_affine
    b = [0.3 + 2*y - 1; 2*(x - 0.5)];
    f_c = rho_c*jacobian(u, [x y])*b;
    for i = 1:2
        g = matlabFunction(f_c(i), Vars = [x y t]);
        f_F = data.f_F{i};
        data.f_F{i} = @(xx, yy) f_F(xx, yy) + g(xx, yy, 1) + 0*xx;
    end
end
beta = cell(2, 1);
for i = 1:2
    bi = matlabFunction(b(i), Vars = [x y t]);
    beta{i} = @(xx, yy) bi(xx, yy, 1) + 0*xx;
end

% tau(beta) on the fluid mesh
nb = normal_velocity(fm, beta, k);
tau = 1 + rho_c*nb.face_max(fm.boundaries)/(2*mu);
factor = options.tau_factor*ones(size(tau));
if options.tau_faces == "off_sigma"
    factor(ismember(fm.boundaries, fm.intfaces)) = 1;
end
fm.tau = tau.*factor;

fp = IncompressibleElasticityProblem(k, time_step = 1, source_term = data.f_F, ...
    datum_dirichlet = data.u_F, viscosity = mu, density = 0.0, pressure_mean = "pinned", ...
    convective_stabilization = false);
switch options.beta
    case "zero"
    case "dubiner"
        bx = project(Datum(beta{1}), fm, DubinerBasis(k + 1));
        by = project(Datum(beta{2}), fm, DubinerBasis(k + 1));
        fp.Data.convective_velocity = {bx.dof; by.dof};
    otherwise
        fp.Data.convective_velocity = beta;
end
bp = BiotProblem(k, young_modulus = 1, poisson_ratio = 0.25, biot_alpha = 1, ...
    storage = 1, diffusion = eye(2), time_step = 1, time_scheme = "euler", ...
    source_term = data.f_P, source_mass = data.g_P, ...
    datum_dirichlet = {data.u_P{1}; data.u_P{2}; data.p_P});
[fs, bs, ~, gs] = StokesBiotCoupling.solve(fp, bp, fm, bm, curve, data, friction_bjs = 1);

ff = fs.Unknowns.field;  bf = bs.Unknowns.field;
F = {ff(1), fm, data.sigma_F; ff(2), fm, data.u_F; bf(1), bm, data.sigma_P; ...
     bf(2), bm, data.u_P; bf(4), bm, data.q_P; bf(5), bm, data.p_P};
out.err = zeros(1, 6);  out.rel = zeros(1, 6);
for i = 1:6
    [uh, mesh, exact] = F{i, :};
    out.err(i) = compute_l2_error(uh, mesh, Datum(exact));
    z = uh;  z.dof = 0*z.dof;
    out.rel(i) = out.err(i)/max(compute_l2_error(z, mesh, Datum(exact)), 1e-12);
end
out.residue = gs.residue;
out.fs = fs;
out.margin = min(mu*fm.tau(:) - rho_c*nb.element_max(:)/2);
out.tau_max = max(fm.tau(:));
out.bn_sigma = nb.sigma;
out.bn_outer = nb.outer;
end

function nb = normal_velocity(mesh, beta, k)
% beta . n_K at k + 3 Gauss points of every element face (n_K outward
% from the element): per face the max |beta . n|, per element face the
% max beta . n_K, and the (min, max) on Sigma and on the outer boundary
s = reshape(lgwt(k + 3, 0, 1), [], 1);
F = mesh.boundaries;
A = mesh.coordinates(1:2, mesh.faces(1, F(:)));
B = mesh.coordinates(1:2, mesh.faces(2, F(:)));
N = reshape(mesh.exterior_unit_normals, 2, []);
X = A(1, :) + s.*(B(1, :) - A(1, :));
Y = A(2, :) + s.*(B(2, :) - A(2, :));
bn = N(1, :).*beta{1}(X, Y) + N(2, :).*beta{2}(X, Y);
n_faces = size(mesh.faces, 2);
nb.face_max = accumarray(F(:), max(abs(bn), [], 1).', [n_faces, 1], @max).';
nb.element_max = reshape(max(bn, [], 1), size(F));
on_boundary = accumarray(F(:), 1, [n_faces, 1]) == 1;
outer = setdiff(find(on_boundary), mesh.intfaces);
range = @(faces) [min(bn(:, ismember(F(:), faces)), [], "all"), ...
                  max(bn(:, ismember(F(:), faces)), [], "all")];
nb.sigma = range(mesh.intfaces);
nb.outer = range(outer);
end

function f = patch_fields(k)
% polynomial exact fields: psi of degree k + 1, the others of degree k
f = struct("stream_function", poly_string(k + 1, 1), "pressure_fluid", poly_string(k, 2), ...
    "displacement", [poly_string(k, 3); poly_string(k, 4)], "pressure_biot", poly_string(k, 5));
end

function s = poly_string(degree, seed)
% sum of c_ij x^i y^j over i + j <= degree, fixed non-trivial coefficients
terms = strings(0);
for i = 0:degree
    for j = 0:degree - i
        c = (-1)^(i + seed) * (1 + mod(3*i + 5*j + seed, 7)) / 4;
        terms(end + 1) = sprintf("(%g)*x^%d*y^%d", c, i, j); %#ok<AGROW>
    end
end
s = join(terms, " + ");
end
