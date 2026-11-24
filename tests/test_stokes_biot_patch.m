% test_stokes_biot_patch.m
% Patch test of the quasi-static Stokes-Biot coupling (StokesBiotCoupling
% .solve) on a flat Sigma with no gap: unit square, Sigma: y = 1/2, fluid
% above, Biot below, fitted_split meshes of sizes (0.8 h, h) with gap 0
% (both meshes fitted to Sigma, Omega_h = Omega).
%
% The exact fields are polynomials that the degree-k spaces contain:
%   stream function of degree k + 1 (u_F of degree k), p_F of degree k,
%   eta and p_P of degree k,
% so sigma_F, sigma_P and q_P have degree k - 1 (constant zero for k = 0).
% StokesBiotCoupling.get_manufactured_solution derives the volume sources
% and the four interface sources from them, so the data satisfy every PDE
% and interface condition. One backward-Euler step from rest (D_t w =
% w/dt), s0 = 1, friction_bjs = 1, Dirichlet data on the whole outer
% boundary of both sides, pinned fluid pressure.
%
% Checked: for k = 1..3 and h = 1/2, 1/4, 1/6, the relative L2 error of
% sigma_F, u_F, sigma_P, eta, q_P, p_P is <= 1e-10 (relative to the L2
% norm of the exact field; absolute where that norm is below 1e-12).
% k = 0 is skipped: the stress-based elements have no curl bubble at
% k = 0, and constant data leave a spurious Darcy flux of order 1e-3.
clc; clear;
tol = 1e-10;
degrees = 1:3;
hs = [1/2, 1/4, 1/6];
names = ["sigma_F", "u_F", "sigma_P", "u_P", "q_P", "p_P"];

n_fail = 0;
for k = degrees
    fields = patch_fields(k);
    for h = hs
        out = stokes_biot_case("flat", k, h, gap = 0, ...
            stream_function = fields.psi, pressure_fluid = fields.p_F, ...
            displacement = fields.eta, pressure_biot = fields.p_P, ...
            storage = 1.0, friction_bjs = 1.0);
        rel = relative_errors(out);
        [worst, at] = max(rel);
        fprintf("k %d  h 1/%-2d  max rel. error %.2e (%s)  residue %.1e\n", ...
            k, round(1/h), worst, names(at), out.residue);
        n_fail = n_fail + check(worst <= tol, ...
            sprintf("k = %d, h = 1/%d: %s relative error %.3e", k, round(1/h), names(at), worst));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_patch:failed", "%d checks failed", n_fail);
end

%=========================================================================%
function f = patch_fields(k)
% polynomial exact fields for degree k: psi of degree k + 1, the others of
% degree k, every monomial present
f.psi = poly_string(k + 1, 1);
f.p_F = poly_string(k, 2);
f.eta = [poly_string(k, 3); poly_string(k, 4)];
f.p_P = poly_string(k, 5);
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

function rel = relative_errors(out)
% L2 errors of the six fields relative to the L2 norms of the exact fields
% (the error of a zero field)
d = out.solutions.data;
ff = out.solutions.fluid.Unknowns.field;
bf = out.solutions.biot.Unknowns.field;
[fm, bm] = out.solutions.meshes{:};
fields = {ff(1), fm, d.sigma_F; ff(2), fm, d.u_F; bf(1), bm, d.sigma_P; ...
          bf(2), bm, d.u_P; bf(4), bm, d.q_P; bf(5), bm, d.p_P};
rel = zeros(1, size(fields, 1));
for i = 1:size(fields, 1)
    [u, mesh, exact] = fields{i, :};
    e = compute_l2_error(u, mesh, Datum(exact));
    z = u;
    z.dof = 0*z.dof;
    scale = compute_l2_error(z, mesh, Datum(exact));
    rel(i) = e / max(scale, 1e-12);
end
end
