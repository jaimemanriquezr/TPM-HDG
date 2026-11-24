% Gallery of TPM mesh pairs: circular and sinusoidal interfaces fitted
% with and without a gap and carved from a background grid, and a flat
% interface with matching and dissimilar structured meshes. Draws the
% meshes, the interface and the transfer paths in one tiled figure.
close all
h = 0.2;

figure();
tiledlayout("flow");

%% Circular interface
r0 = .25;
x0 = [.5; .5];
circle = CircularCurve(x0, r0);
draw_mesh_couplings(h, circle, @(I, h) fitted_split(I, h, 0), nexttile)
proper_fit_closed = @(I, h) fitted_split(I, h, 0.1);
draw_mesh_couplings(h, circle, proper_fit_closed, nexttile)
% Cut elements discarded: both subdomains stay strictly inside their exact
% domains and the removed band is the gap. Both options are the defaults,
% named here for illustration. The crossed background grid admits more
% face directions, so the pruned staircase can meet the 45 degree
% feasibility bound.
cut_discarded = @(I, h) StructuredSimplicialGrid(1/h, 1/h, "crossed").split(I, ...
    cut_elements="discard", max_normal_deviation=45);
draw_mesh_couplings(h, circle, cut_discarded, nexttile)

%% Sinusoidal interface
sine = SmoothCurveHorizontal(".5 + .1*sin(2*pi*t)");
draw_mesh_couplings(h, sine, @(I, h) fitted_split(I, h, 0), nexttile)
draw_mesh_couplings(h, sine, @(I, h) fitted_split(I, h, .1), nexttile)
draw_mesh_couplings(h, sine, cut_discarded, nexttile)

%% Flat interface
y0 = FlatCurveHorizontal(0.5);
flat_split = @(N) @(I, h) fitted_split(I, [h, N*h], .1, type="structured");
draw_mesh_couplings(h, y0, flat_split(1), nexttile);
draw_mesh_couplings(h, y0, flat_split(2), nexttile);
