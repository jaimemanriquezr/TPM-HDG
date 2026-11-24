%-------------------------------------------------------------------------
%Description: InterfaceProblem.stokes_stokes (L-u-p formulation) with a
%             pressure whose normal derivative on Sigma does not vanish.
%
%             Sigma is the line y = y0 = 1/3, fitted meshes retreating
%             h^2 (fitted_split), both sides with the constrained
%             pressure closure. The manufactured pressure
%             p = cos(pi x) sin(pi y) has zero x-average, so the
%             zero-mean closure is exact on any horizontal strip, and
%             d_n p = pi cos(pi x) cos(pi y0) on Sigma, nonzero. Two cases:
%               "equal nu":     nu_1 = nu_2 = 0.7, psi = sin(pi x) sin(pi y),
%                               so E(L_i) n_i enters T1 nontrivially;
%               "different nu": nu_1 = 0.7 (upper), nu_2 = 1.3 (lower),
%                               psi = sin(pi x) (1 + (y - y0)^3), so
%                               d_y u = 0 on Sigma and the traction
%                               nu_i (grad u) n - p n is continuous for
%                               any nu_1, nu_2 (checked symbolically).
%             u = (-psi_y, psi_x), f_i = -nu_i lap u + grad p per side.
%             Checks: mean EOC of u, L = grad u and p >= k + 0.5 on
%             h = 1/4, 1/8, 1/16, k = 1, 2.
%             Level-free side: "different nu" data, k = 1, with
%             p = cos(pi x) sin(pi y) + y (strip means differ), the upper
%             side constrained to zero mean and the lower side without
%             the closure; its level is fixed only through T1. Checks:
%             mean EOC of the lower p >= k + 0.5 against p - mean_upper(p)
%             and its mean error decreasing below 1e-4.
%-------------------------------------------------------------------------
n_failed = 0;

y0 = 1/3;
interface = FlatCurveHorizontal(y0);
syms x y

cases = struct( ...
    "name", {"equal nu", "different nu"}, ...
    "nu", {[0.7, 0.7], [0.7, 1.3]}, ...
    "psi", {sin(pi*x) .* sin(pi*y), sin(pi*x) .* (1 + (y - y0)^3)});
p_sym = cos(pi*x) .* sin(pi*y);
dnp = simplify(subs(diff(p_sym, y), y, y0));
fprintf("d_n p on Sigma = %s, max |d_n p| = %.4f\n", string(dnp), ...
    double(max(abs(subs(dnp, x, linspace(0, 1, 101))))));

for c = 1:numel(cases)
    nu = cases(c).nu;
    psi = cases(c).psi;
    u_sym = [-diff(psi, y); diff(psi, x)];
    L_sym = jacobian(u_sym, [x y]);
    % traction balance across Sigma (n = e_y on both sides' common line)
    jump = simplify(subs((nu(1) - nu(2)) * L_sym * [0; 1], y, y0));
    fprintf("\n=== %s (nu_1 = %.1f, nu_2 = %.1f): traction jump on Sigma = [%s; %s]\n", ...
        cases(c).name, nu(1), nu(2), string(jump(1)), string(jump(2)));
    n_failed = n_failed + check(cases(c).name + ": traction continuous", all(jump == 0));

    u_ = Datum(u_sym);
    L_ = Datum(L_sym);
    p_ = Datum(p_sym);
    f_han = cell(1, 2);
    for j = 1:2
        f_sym = -nu(j)*[laplacian(u_sym(1), [x y]); laplacian(u_sym(2), [x y])] ...
              + gradient(p_sym, [x y]);
        f_han{j} = {matlabFunction(f_sym(1), "Vars", [x y]); ...
                    matlabFunction(f_sym(2), "Vars", [x y])};
    end

    for k = [1, 2]
        Pk = DubinerBasis(k);
        n_ref = 3;
        [err_u, err_L, err_p, h_list] = deal(zeros(n_ref, 1));
        for i = 1:n_ref
            h = 1/4 * (1/2)^(i-1);
            h_list(i) = h;
            [upper_mesh, lower_mesh] = fitted_split(interface, h, h^2);
            ms = {upper_mesh, lower_mesh};
            pb = cell(1, 2);
            for j = 1:2
                pb{j} = FEProblem( ...
                    Spaces = [FESpace([Pk.to_tensor(), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), ...
                              FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])], ...
                    Order = k);
                pb{j}.Mesh = ms{j};
                [lhs, flux, rhs] = StokesMatrices(k, ms{j}, nu(j), f_han{j});
                pb{j}.Operator = AffineOperator(lhs + rhs);
                pb{j}.Operator = pb{j}.Operator + (flux);
                pb{j}.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");
            end

            coupled = InterfaceProblem.stokes_stokes(pb{1}, pb{2}, interface, nu(1), nu(2));
            [sol1, sol2, gs] = solve(coupled, k);

            sols = {sol1, sol2};
            [eu, eL, ep] = deal(zeros(1, 2));
            for j = 1:2
                eL(j) = compute_l2_error(sols{j}.Unknowns.field(1), ms{j}, L_);
                eu(j) = compute_l2_error(sols{j}.Unknowns.field(2), ms{j}, u_);
                ep(j) = compute_l2_error(sols{j}.Unknowns.field(3), ms{j}, p_);
            end
            err_u(i) = max(eu);  err_L(i) = max(eL);  err_p(i) = max(ep);
            fprintf("k = %d, h = %.4f: err u = (%.3e, %.3e), L = (%.3e, %.3e), p = (%.3e, %.3e), residue = %.1e\n", ...
                k, h, eu, eL, ep, gs.residue);
        end
        rate = @(e) compute_eoc(e, h_list);
        eoc = struct("u", rate(err_u), "L", rate(err_L), "p", rate(err_p));
        for fld = ["u", "L", "p"]
            fprintf("k = %d: EOC(%s) = %s (expected ~ %d)\n", ...
                k, fld, num2str(eoc.(fld).', "%.2f "), k + 1);
            n_failed = n_failed + check(sprintf("%s, k = %d: mean EOC(%s) >= k + 0.5", ...
                cases(c).name, k, fld), mean(eoc.(fld)) >= k + 0.5);
        end
    end
end

% one side constrained (upper), the other level-free (lower)
nu = cases(2).nu;
u_sym = [-diff(cases(2).psi, y); diff(cases(2).psi, x)];
u_ = Datum(u_sym);
p_lf = p_sym + y;
f_han = cell(1, 2);
for j = 1:2
    f_sym = -nu(j)*[laplacian(u_sym(1), [x y]); laplacian(u_sym(2), [x y])] ...
          + gradient(p_lf, [x y]);
    f_han{j} = {matlabFunction(f_sym(1), "Vars", [x y]); ...
                matlabFunction(f_sym(2), "Vars", [x y])};
end
k = 1;
Pk = DubinerBasis(k);
n_ref = 3;
[err_p, err_mean, h_list] = deal(zeros(n_ref, 1));
% mean of p_lf over a [0,1] x [a, b] mesh is (a + b)/2
mean_of = @(mesh) (min(mesh.vertices(2, :, :), [], "all") + max(mesh.vertices(2, :, :), [], "all"))/2;
fprintf("\n=== level-free lower side (nu_1 = %.1f, nu_2 = %.1f), p = %s\n", nu, string(p_lf));
for i = 1:n_ref
    h = 1/4 * (1/2)^(i-1);
    h_list(i) = h;
    [upper_mesh, lower_mesh] = fitted_split(interface, h, h^2);
    ms = {upper_mesh, lower_mesh};
    pb = cell(1, 2);
    for j = 1:2
        global_spaces = [LegendreBasis(k).to_vector(2), DubinerBasis(0)];
        [lhs, flux, rhs] = StokesMatrices(k, ms{j}, nu(j), f_han{j});
        if j == 1
            global_spaces = [global_spaces, ConstantBasis()]; %#ok<AGROW>
        else
            flux.Arrays = flux.Arrays(~cellfun(@(A) any(A.id == 7), flux.Arrays));
        end
        pb{j} = FEProblem( ...
            Spaces = [FESpace([Pk.to_tensor(), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), ...
                      FESpace(global_spaces)], ...
            Order = k);
        pb{j}.Mesh = ms{j};
        pb{j}.Operator = AffineOperator(lhs + rhs);
        pb{j}.Operator = pb{j}.Operator + (flux);
        pb{j}.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");
    end
    [~, sol2, gs] = solve(InterfaceProblem.stokes_stokes(pb{1}, pb{2}, interface, nu(1), nu(2)), k);

    p_ref = p_lf - mean_of(upper_mesh);
    area = sum(lower_mesh.volume);
    e0 = compute_l2_error(sol2.Unknowns.field(3), lower_mesh, Datum(p_ref));
    e1 = compute_l2_error(sol2.Unknowns.field(3), lower_mesh, Datum(p_ref + 1));
    err_p(i) = e0;
    err_mean(i) = abs(e0^2 + area - e1^2)/(2*area);   % |mean(p_h - p_ref)|
    fprintf("h = %.4f: lower p error = %.3e, exact lower mean = %+.4f, mean error = %.2e, residue = %.1e\n", ...
        h, e0, mean_of(lower_mesh) - mean_of(upper_mesh), err_mean(i), gs.residue);
end
eoc_p = compute_eoc(err_p, h_list);
fprintf("k = %d: EOC(p, level-free side) = %s (expected ~ %d)\n", k, num2str(eoc_p.', "%.2f "), k + 1);
n_failed = n_failed + check("level-free side: mean EOC(p) >= k + 0.5", mean(eoc_p) >= k + 0.5);
n_failed = n_failed + check("level-free side: mean error decreasing, < 1e-4", ...
    all(diff(err_mean) < 0) && err_mean(end) < 1e-4);

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
