%TEST_STOKES_BIOT_PHYSICAL_TAU  Patch test of Stokes-Biot with a solid
% stabilization tau_s = stabilization_solid * mesh.tau far from 1.
%
% E = 3e5 and tau_s = 2 mu_P (the physical rule with L = 1), flat Sigma,
% gap h^2, k = 1..3, h = 1/2, 1/4: the polynomial patch solution of
% test_stokes_biot_patch is reproduced in every field to a relative L2
% error <= 1e-8. The tolerance is not the 1e-10 of
% test_stokes_biot_patch: with E = 3e5 the k = 3 round-off reaches about
% 2e-9 whatever tau_s is (the same with tau_s = 1 and 1e9).
%
% What it does not check: at the exact solution eta = eta_hat on the
% faces, so the tau_s terms vanish and a mismatch of tau_s between the
% element (BiotMatrices) and the rows (stokes_biot_rows) passes too.

degrees = 1:3;
hs = [1/2, 1/4];
tol = 1e-8;
E = 3e5;  nu = 0.25;
[~, mu_p] = LinearElasticityProblem.get_lame_parameters(E, nu);
names = ["sigma_F", "u_F", "sigma_P", "u_P", "q_P", "p_P"];

n_fail = 0;
for k = degrees
    f = patch_fields(k);
    for h = hs
        out = stokes_biot_case("flat", k, h, ...
            stream_function = f.psi, pressure_fluid = f.p_F, ...
            displacement = f.eta, pressure_biot = f.p_P, ...
            young_modulus = E, poisson_ratio = nu, ...
            stabilization_solid = 2*mu_p, storage = 1.0, friction_bjs = 1.0);
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
    error("test_stokes_biot_physical_tau:failed", "%d checks failed", n_fail);
end

%=========================================================================%
function f = patch_fields(k)
% the polynomial fields of test_stokes_biot_patch
f.psi = poly_string(k + 1, 1);
f.p_F = poly_string(k, 2);
f.eta = [poly_string(k, 3); poly_string(k, 4)];
f.p_P = poly_string(k, 5);
end

function s = poly_string(degree, seed)
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
