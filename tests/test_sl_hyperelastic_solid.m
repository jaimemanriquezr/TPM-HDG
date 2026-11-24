%TEST_SL_HYPERELASTIC_SOLID  Patch tests of the semi-Lagrangian Eulerian hyperelastic solid.
%
%   The transport alone: the t = 0 patch level carried once by the affine
%   (rotation + translation) feet of v = V0 + omega (x - c0)^perp matches
%   the projections of the exactly composed fields, on a fixed and on a
%   translated carve. The march: a uniform material velocity and u = U(x -
%   V0 t, t), U of degree k (linear law) or 1 (neo-Hookean), quadratic in
%   t, are reproduced on a fixed and on a translating carve by BDF2 and by
%   CN with the transported residual, k = 1, 2, 3.
%
%See also: sl_hyperelastic_solid_march
failures = 0;
tol = 1e-10;
for k = 1:3
    for A = [0 0.1]
        o = sl_hyperelastic_solid_march(k, 8, A, 0.05, 0, fields = "patch", transport_check = 2);
        fprintf("transport k = %d, A = %.1f: relative error %.3e (uncovered %.3f)\n", k, A, ...
            o.transport_error, o.uncovered);
        failures = failures + ~(o.transport_error <= tol);
    end
end
cases = {"linear", [], "bdf2"; "linear", [], "cn"; "neo_hookean", 1, "bdf2"; "neo_hookean", 1, "cn"};
for k = 1:3
    for c = 1:size(cases, 1)
        for A = [0 0.1]
            for N = [8 16]
                o = sl_hyperelastic_solid_march(k, N, A, 0.05, 3, fields = "patch", law = cases{c, 1}, ...
                    patch_degree = cases{c, 2}, scheme = cases{c, 3});
                e = max([o.e_v/o.n_v, o.e_L/o.n_L, o.e_sigma/o.n_sigma, o.e_u/o.n_u]);
                fprintf("patch %s %s k = %d, A = %.1f, N = %d: relative error %.3e, Newton %s\n", ...
                    cases{c, 1}, cases{c, 3}, k, A, N, e, mat2str(o.iterations.'));
                failures = failures + ~(e <= tol && all(o.converged));
            end
        end
    end
end
if failures == 0
    fprintf("test_sl_hyperelastic_solid: all checks passed\n");
else
    fprintf("test_sl_hyperelastic_solid: %d check(s) failed\n", failures);
end
