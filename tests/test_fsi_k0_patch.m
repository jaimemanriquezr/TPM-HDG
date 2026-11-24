% test_fsi_k0_patch.m
% Degree-0 patch test of the FSI coupling (FluidStructureCoupling.solve,
% steady, time = Inf). Constant data that P0 reproduces exactly: fluid
% stream function x (u_F = (0, 1)), zero pressure, solid displacement
% u_S = (1, 0), so every stress vanishes. Checked, at machine precision on
% sigma, u, gamma of both sides:
%   - flat interface, fitted meshes, densities (1, 1), pinned and
%     constrained pressure closures;
%   - flat interface, fitted meshes, densities (0, 0);
%   - sine interface, interior (unfitted, cut elements discarded) meshes,
%     densities (1, 1).
% With DubinerBasis(0)/LegendreBasis(0) the face matrices are 1-by-1 per
% face and PEERSBasis(0) has no bubble modes; both shapes are exercised.
clc; clear;
k = 0; h = 1/8; tol = 1e-11;
flat = FlatCurveHorizontal(0.5);
sine = SmoothCurveHorizontal("0.5 + .05*sin(2*pi*t)", [0 1]);
[flat_fluid, flat_solid] = fitted_split(flat, [h h], 0);
[sine_fluid, sine_solid] = interface_meshes(sine, "interior", h, max_normal_deviation = 45);

cases = {
    "flat, fitted, rho = (1, 1), pinned",      flat, flat_fluid, flat_solid, 1, "pinned"
    "flat, fitted, rho = (1, 1), constrained", flat, flat_fluid, flat_solid, 1, "constrained"
    "flat, fitted, rho = (0, 0), pinned",      flat, flat_fluid, flat_solid, 0, "pinned"
    "sine, interior, rho = (1, 1), pinned",    sine, sine_fluid, sine_solid, 1, "pinned"};

n_fail = 0;
for c = 1:size(cases, 1)
    [what, curve, fluid_mesh, solid_mesh, rho, closure] = cases{c, :};
    [fp, sp, sources, pd] = FluidStructureCoupling.load_problems(k, curve, Inf, 1.0, ...
        "x", "0", ["1"; "0"], rho, 1.0, rho, 1.0, 0.25, 0.0, pressure_mean = closure);
    [fs, ss] = FluidStructureCoupling.solve(fp, sp, fluid_mesh, solid_mesh, curve, ...
        sources{:}, "steady", patch_weight = "jacobian");
    eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
    eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
    e = max([eF(:); eS(:)]);
    fprintf("%-42s max error %.3e\n", what, e);
    n_fail = n_fail + check(e < tol, sprintf("%s: max error %.3e", what, e));
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_k0_patch:failed", "%d checks failed", n_fail);
end
