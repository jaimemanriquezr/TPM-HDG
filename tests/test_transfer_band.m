% test_transfer_band.m
% Verification of transfer_band_jacobian, the signed area element of the band
% between an immersed mesh's boundary and the true boundary.
%
% The band integral of 1 is an exact analytic quantity:
%
%     int_band 1 dA  ==  area(true domain) - sum_K |K|
%
% immerse_mesh discards every cut element, so the computational polygon lies
% strictly inside the disk and the deficit is positive.
%
% The identity holds up to the s-quadrature error (xbar(s) is
% transcendental for a circle); section 2 shows that error reaching roundoff
% with a handful of nodes. Section 3 compares finite-difference and analytic
% mapping Jacobians.
clc; clear;

n_fail = 0;

disk = @(c, R) struct( ...
    "phi",     @(X) vecnorm(X - c, 2, 1) - R, ...
    "mapping", @(X) c + R * (X - c) ./ vecnorm(X - c, 2, 1), ...
    "area",    pi * R^2, ...
    "c", c, "R", R);

band_area = @(J, info) sum(J .* info.weights, 'all');

cases = { ...
    disk([0; 0],       0.75), "centred R = 0.75"; ...
    disk([0; 0],       0.50), "centred R = 0.50"; ...
    disk([0.11; -0.07], 0.60), "off-centre R = 0.60"};

%% ------------------------------------------------------------------ %%
%% 1. area deficit under refinement
%% ------------------------------------------------------------------ %%
fprintf("1. band integral of 1 vs the analytic area deficit (8 s-nodes)\n\n");

for ci = 1:size(cases, 1)
    g = cases{ci, 1};
    fprintf("   %s\n", cases{ci, 2});

    for n = [8 16 32 64]
        mesh = immerse_mesh(StructuredSimplicialGrid(n, n, -1, 1, -1, 1), g.phi);

        [J, info] = transfer_band_jacobian(mesh, g.mapping, ...
            quadrature_face = 8, quadrature_path = 2);

        deficit = g.area - sum(mesh.volume);
        computed = band_area(J, info);
        rel = abs(computed - deficit) / abs(deficit);

        fprintf("     n = %3d  deficit = %+.10e  band = %+.10e  rel = %.2e  under-covers = %d\n", ...
            n, deficit, computed, rel, info.under_covers);

        if rel > 1e-10
            fprintf("       FAILED: band area does not match the analytic deficit\n");
            n_fail = n_fail + 1;
        end
        if ~info.under_covers
            fprintf("       FAILED: paths should run outward (polygon strictly inside)\n");
            n_fail = n_fail + 1;
        end
        if deficit <= 0
            fprintf("       FAILED: deficit should be positive for a discarded-cut mesh\n");
            n_fail = n_fail + 1;
        end
    end
    fprintf("\n");
end

%% ------------------------------------------------------------------ %%
%% 2. the s-quadrature floor
%% ------------------------------------------------------------------ %%
% t enters det DX linearly, so 2 path nodes are exact in t for any geometry;
% only the s-rule is approximate. The error should fall steeply with s-nodes
% and then sit on roundoff.
fprintf("2. dependence on the number of s-nodes (n = 32, centred R = 0.75)\n\n");

g = cases{1, 1};
mesh = immerse_mesh(StructuredSimplicialGrid(32, 32, -1, 1, -1, 1), g.phi);
deficit = g.area - sum(mesh.volume);

% Run it both ways: the finite-difference xbar' has its own accuracy floor
% (about eps^(2/3) for a central difference), which the analytic column
% separates from the quadrature error.
analytic = @(X) radial_jacobian(X, g.c, g.R);
fprintf("     %-9s %-14s %-14s\n", "s-nodes", "FD xbar'", "analytic xbar'");
for ns = [1 2 3 5 8 12]
    [Jf, inf_f] = transfer_band_jacobian(mesh, g.mapping, ...
        quadrature_face = ns, quadrature_path = 2);
    [Ja, inf_a] = transfer_band_jacobian(mesh, g.mapping, ...
        quadrature_face = ns, quadrature_path = 2, mapping_jacobian = analytic);
    rel_f = abs(band_area(Jf, inf_f) - deficit) / abs(deficit);
    rel_a = abs(band_area(Ja, inf_a) - deficit) / abs(deficit);
    fprintf("     %-9d %-14.3e %-14.3e\n", ns, rel_f, rel_a);
end

% The quadrature is resolved by a handful of nodes; check the analytic route
% actually reaches roundoff, so the FD plateau above is attributable.
[Ja, inf_a] = transfer_band_jacobian(mesh, g.mapping, ...
    quadrature_face = 12, quadrature_path = 2, mapping_jacobian = analytic);
rel_a = abs(band_area(Ja, inf_a) - deficit) / abs(deficit);
if rel_a > 1e-14
    fprintf("       FAILED: with an exact xbar' the band area should reach " + ...
        "roundoff (got %.2e)\n", rel_a);
    n_fail = n_fail + 1;
end
fprintf("\n");

% and that the path rule really is exact at 2 nodes
[J2, i2] = transfer_band_jacobian(mesh, g.mapping, ...
    quadrature_face = 8, quadrature_path = 2);
[J9, i9] = transfer_band_jacobian(mesh, g.mapping, ...
    quadrature_face = 8, quadrature_path = 9);
drift = abs(band_area(J9, i9) - band_area(J2, i2)) / abs(deficit);
fprintf("     path nodes 2 vs 9: rel difference = %.2e  (det is linear in t)\n\n", drift);
if drift > 1e-13
    fprintf("       FAILED: the path rule is not exact at 2 nodes\n");
    n_fail = n_fail + 1;
end

%% ------------------------------------------------------------------ %%
%% 3. finite-difference xbar' against the analytic map Jacobian
%% ------------------------------------------------------------------ %%
% For xbar = c + R (x-c)/r, D(xbar) = (R/r)(I - uu'), u = (x-c)/r. The band
% area computed either way must agree far below the quadrature error above.
fprintf("3. finite-difference vs analytic mapping Jacobian\n\n");

for ci = 1:size(cases, 1)
    g = cases{ci, 1};
    mesh = immerse_mesh(StructuredSimplicialGrid(32, 32, -1, 1, -1, 1), g.phi);

    analytic = @(X) radial_jacobian(X, g.c, g.R);

    [J_fd, i_fd] = transfer_band_jacobian(mesh, g.mapping, ...
        quadrature_face = 8, quadrature_path = 2);
    [J_an, i_an] = transfer_band_jacobian(mesh, g.mapping, ...
        quadrature_face = 8, quadrature_path = 2, mapping_jacobian = analytic);

    rel = abs(band_area(J_an, i_an) - band_area(J_fd, i_fd)) ...
        / abs(band_area(J_an, i_an));
    node_diff = max(abs(J_an - J_fd), [], 'all') / max(abs(J_an), [], 'all');

    fprintf("   %-22s  area rel diff = %.2e   worst node rel diff = %.2e\n", ...
        cases{ci, 2}, rel, node_diff);

    if rel > 1e-9
        fprintf("       FAILED: FD and analytic xbar' disagree on the band area\n");
        n_fail = n_fail + 1;
    end
end

%% ------------------------------------------------------------------ %%
%% 4. over-covering and folded bands: the integral stays signed
%% ------------------------------------------------------------------ %%
% Carving a larger disk keeps elements beyond the true circle (negative
% deficit); carving a shifted or slightly larger disk makes the polygon
% cross the circle, so the ruled surface folds and J takes both signs. In
% every case the band integral is the signed deficit (Green's theorem).
fprintf("\n4. over-covering and folded meshes (n = 32, centred R = 0.60)\n\n");

R = 0.60; n = 32; h = 2/n;
g = disk([0; 0], R);
background = StructuredSimplicialGrid(n, n, -1, 1, -1, 1);
carves = { ...
    @(X) vecnorm(X, 2, 1) - (R + 3*h),        "over-covering (carve R + 3h)",   false; ...
    @(X) vecnorm(X - [0.05; 0], 2, 1) - R,    "folded (carve shifted by 0.05)", true; ...
    @(X) vecnorm(X, 2, 1) - (R + 0.7*h),      "folded (carve R + 0.7h)",        true};
for ci = 1:size(carves, 1)
    mesh = immerse_mesh(background, carves{ci, 1});
    [J, info] = transfer_band_jacobian(mesh, g.mapping, ...
        quadrature_face = 12, quadrature_path = 2);

    deficit = g.area - sum(mesh.volume);
    computed = band_area(J, info);
    rel = abs(computed - deficit) / abs(deficit);
    folded = any(J < 0, 'all') && any(J > 0, 'all');

    fprintf("   %-32s deficit = %+.10e  band = %+.10e  rel = %.2e  folded = %d\n", ...
        carves{ci, 2}, deficit, computed, rel, folded);

    if rel > 1e-10
        fprintf("       FAILED: band area does not match the signed deficit\n");
        n_fail = n_fail + 1;
    end
    if folded ~= carves{ci, 3}
        fprintf("       FAILED: expected folded = %d\n", carves{ci, 3});
        n_fail = n_fail + 1;
    end
end
mesh = immerse_mesh(background, carves{1, 1});
if g.area - sum(mesh.volume) >= 0
    fprintf("       FAILED: the over-covering deficit should be negative\n");
    n_fail = n_fail + 1;
end

fprintf("\n%d check(s) failed.\n", n_fail);

% ----------------------------------------------------------------------- %
function D = radial_jacobian(X, c, R)
%RADIAL_JACOBIAN  D of x -> c + R (x-c)/|x-c|, as (2 x 2 x N).
V = X - c;
r = vecnorm(V, 2, 1);
u = reshape(V ./ r, 2, 1, []);
D = (R ./ reshape(r, 1, 1, [])) .* (eye(2) - pagemtimes(u, 'none', u, 'transpose'));
end
