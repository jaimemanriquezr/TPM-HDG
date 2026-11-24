% test_stokes_biot_patch_gap.m
% Stokes-Biot (HDG + TPM): coupled patch test on the flat interface
% Sigma: y = 1/2 with a gap between the meshes (fluid above, Biot below,
% fitted_split sizes (0.8 h, h) through stokes_biot_case).
%
% Every exact field is a polynomial of degree <= k: the stream function
% has degree k + 1 (so u_F has degree k), and p_F, eta and p_P have degree
% k. The fluid and Biot discrete spaces contain them, and the Taylor
% transfer along the paths is exact for polynomials of degree <= k, so the
% coupled solution must reproduce them up to round-off: every L2 error
% (sigma_F, u_F, sigma_P, u_P, q_P, p_P) relative to the L2 norm of the
% discrete field <= max(1e-10, 100 r) (absolute when that norm is below
% 1e-12), r the round-off of the global solve measured by one step of
% iterative refinement (solve_roundoff). With the gap h the system's
% condition number reaches 1e9 to 1e10 at k = 3, and the error, about 4 r,
% exceeds 1e-10 (1.3e-10 at h = 1/4) while staying at round-off. One
% Euler step from rest (Delta t = 1), s0 = 1, friction_bjs = 1,
% all-Dirichlet outer boundaries; the interface sources
% of StokesBiotCoupling.get_manufactured_solution carry the jumps of the
% independent fields.
%
%   gaps  delta = h^2 and delta = h
%   h     1/4, 1/8, 1/16
%   k     1, 2, 3 (k = 0 is skipped: the stress-based elements have no
%         curl bubble at k = 0, and the rotation, and in the fluid the
%         stress and the pressure, do not converge there)
clc; clear;

tol = 1e-10;
degrees = 1:3;
hs = [1/4, 1/8, 1/16];
gaps = ["h2", "h"];
fields = ["sigma_F", "u_F", "sigma_P", "u_P", "q_P", "p_P"];
% the discrete field of out.solutions behind each error: (side, index)
side = ["fluid", "fluid", "biot", "biot", "biot", "biot"];
index = [1, 2, 1, 2, 4, 5];
base = {'storage', 1.0, 'friction_bjs', 1.0};

n_fail = 0;
worst = zeros(numel(gaps), numel(degrees));
for g = 1:numel(gaps)
    for kk = 1:numel(degrees)
        k = degrees(kk);
        poly = polynomial_fields(k);
        for h = hs
            delta = h^2;
            if gaps(g) == "h", delta = h; end
            out = stokes_biot_case("flat", k, h, base{:}, poly{:}, 'gap', delta);
            rel = zeros(1, numel(fields));
            for j = 1:numel(fields)
                u_h = out.solutions.(side(j)).Unknowns.field(index(j));
                mesh = out.solutions.meshes{1 + (side(j) == "biot")};
                scale = compute_l2_error(u_h, mesh);
                if scale < 1e-12, scale = 1; end
                rel(j) = out.errors.(fields(j)) / scale;
            end
            worst(g, kk) = max(worst(g, kk), max(rel));
            roundoff = solve_roundoff(out.solutions.global_system);
            bound = max(tol, 100*roundoff);
            fprintf("gap %-2s k %d h 1/%-2d residue %.1e  solve round-off %.1e  %s\n", gaps(g), k, ...
                round(1/h), out.residue, roundoff, strjoin(compose("%s %.1e", fields.', rel.'), "  "));
            for j = 1:numel(fields)
                n_fail = n_fail + check(rel(j) <= bound, sprintf( ...
                    "gap %s, k = %d, h = 1/%d: relative error of %s %.2e > %.1e", ...
                    gaps(g), k, round(1/h), fields(j), rel(j), bound));
            end
        end
    end
end
for g = 1:numel(gaps)
    fprintf("gap %-2s worst relative error per k = %s: %s\n", gaps(g), ...
        mat2str(degrees), strjoin(compose("%.1e", worst(g, :)), ", "));
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_patch_gap:failed", "%d checks failed", n_fail);
end

function options = polynomial_fields(k)
% exact fields of degree <= k (stream function k + 1), as the name-value
% options of stokes_biot_case
options = {'stream_function', sprintf("(0.3*x - 0.7*y + 0.2)^%d + 0.5*(0.6*x + 0.4*y - 0.1)^%d", k + 1, k + 1), ...
    'pressure_fluid', sprintf("0.4*(0.4*x + 0.6*y - 0.1)^%d + 0.3", k), ...
    'displacement', [sprintf("(0.5*x + 0.2*y + 0.1)^%d", k); sprintf("0.2 - (0.1*x - 0.4*y + 0.3)^%d", k)], ...
    'pressure_biot', sprintf("(0.2*x + 0.5*y + 0.4)^%d - 0.1", k)};
end
