% test_stokes_biot_patch_circle.m
% Patch test of the quasi-static Stokes-Biot coupling on a curved
% interface (StokesBiotCoupling.solve through stokes_biot_case): Sigma the
% circle of centre (1/2, 1/2) and radius 1/4, Biot inside, fluid outside,
% interface_meshes(curve, "fitted", 0.8 h, gap = h^2, h_solid = h), so
% neither mesh reaches Sigma and both are polygons.
%
% For k = 1..3 every exact field is a full polynomial of degree <= k
% (k = 0 is skipped: the stress elements carry no curl bubble at k = 0,
% and the rotations, fluid stress and pressure do not converge there):
% stream function of degree k + 1 (u_F of degree k), p_F, eta and p_P of
% degree k, with fixed non-zero coefficients; hence sigma_F, sigma_P, q_P
% are of degree <= k too. The interface sources T_1..T_4 of
% StokesBiotCoupling.get_manufactured_solution carry the exact normal of
% the circle (not a polynomial) and are evaluated at the patch points on
% Sigma. Storage s0 = 1, friction_bjs = 1, one Euler step from rest.
% The exact solution lies in the discrete space and the Taylor transfer
% along the paths is exact for it, so the coupled solution must reproduce
% it: relative L2 error <= 1e-10 for sigma_F, u_F, sigma_P, eta, q_P,
% p_P (relative to the L2 norm of the discrete field; a field whose exact
% value vanishes is measured against the largest field norm) on h = 1/8,
% 1/12, 1/16.
clc; clear;
addpath(fullfile(fileparts(mfilename("fullpath")), "helpers"));
tol = 1e-10;
hs = [1/8, 1/12, 1/16];
names = ["sigma_F", "u_F", "sigma_P", "u_P", "q_P", "p_P"];

n_fail = 0;
for k = 1:3
    psi = polynomial(k + 1, 1);
    p_F = polynomial(k, 2);
    eta = [polynomial(k, 3); polynomial(k, 4)];
    p_P = polynomial(k, 5);
    for h = hs
        out = stokes_biot_case("circle", k, h, stream_function = psi, ...
            pressure_fluid = p_F, displacement = eta, pressure_biot = p_P, ...
            storage = 1.0, friction_bjs = 1.0, gap = h^2);
        [fluid_mesh, biot_mesh] = out.solutions.meshes{:};
        ff = out.solutions.fluid.Unknowns.field;
        bf = out.solutions.biot.Unknowns.field;
        norms = [compute_l2_error(ff(1), fluid_mesh), compute_l2_error(ff(2), fluid_mesh), ...
            compute_l2_error(bf(1), biot_mesh), compute_l2_error(bf(2), biot_mesh), ...
            compute_l2_error(bf(4), biot_mesh), compute_l2_error(bf(5), biot_mesh)];
        errors = arrayfun(@(name) out.errors.(name), names);
        % a field whose exact value vanishes is measured against the
        % largest field norm
        scale = norms;
        scale(norms < 1e-8*max(norms)) = max(norms);
        relative = errors ./ scale;
        [worst, j] = max(relative);
        fprintf("k %d  h 1/%-3d  max relative error %.2e (%s)  residue %.1e\n", ...
            k, round(1/h), worst, names(j), out.residue);
        n_fail = n_fail + check(worst <= tol, sprintf( ...
            "k %d, h 1/%d: relative %s error %.3e", k, round(1/h), names(j), worst));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_patch_circle:failed", "%d checks failed", n_fail);
end

function s = polynomial(n, seed)
% a full polynomial of degree n in x, y with fixed non-zero coefficients
terms = string(seed/7);
for m = 1:n
    for i = 0:m
        c = (1 + mod(seed*(i + 3) + m, 5))/(2 + i + m);
        if mod(i + seed, 2), c = -c; end
        terms(end + 1) = sprintf("%.17g*x^%d*y^%d", c, m - i, i); %#ok<AGROW>
    end
end
s = join(terms, " + ");
end
