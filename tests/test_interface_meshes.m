% interface_meshes invariants per family: n_shared equals all interface
% faces for the protrude families (both sides are the grid), 0 for
% interior, and at least 90% of them for fitted (mesh2d may refine a thin
% polygon's boundary on one side); area_gap = 0 except for interior.
% interface_polygon on the fitted mesh must close on a closed curve and
% stay open on a graph.
clc; clear;
n_fail = 0;
cases = {"kidney",  discretize(KidneyCurve(), n_grid = 4096),  [-1 1 -1 1],     "outer";
         "airfoil", AirfoilCurve(-1/10, 1/10, 0.02),           [-3 3 -0.4 0.6], "outer";
         "circle",  CircularCurve([.5; .5], .4),               [0 1 0 1],       "inner";
         "curved",  SmoothCurveHorizontal("0.5 + .05*sin(2*pi*t)"), [0 1 0 1],  "lower"};
h = 1/12;
for c = 1:size(cases, 1)
    [name, curve, box, side] = cases{c, :};
    fm_fitted = [];
    for family = ["fitted", "interior", "protrude_fluid", "protrude_solid"]
        if name == "curved" && family ~= "fitted", continue; end   % graph curves: fitted only here
        [fm, sm, info] = interface_meshes(curve, family, h, box = box, fluid_side = side);
        nF = numel(fm.intfaces);  nS = numel(sm.intfaces);
        % fitted pairs are zero-gap but mesh2d may refine a thin polygon's
        % boundary on one side (the airfoil), so conformity is asserted
        % only for the split families, where both sides are the grid
        ok_shared = (family == "interior" && info.n_shared == 0) || ...
                    (family == "fitted"   && info.n_shared >= 0.9*min(nF, nS)) || ...
                    (ismember(family, ["protrude_fluid", "protrude_solid"]) ...
                        && info.n_shared == nF && nF == nS);
        ok_gap = (family == "interior") == (abs(info.area_gap) > 1e-10);
        fprintf("%-8s %-15s el %5d/%5d  ifaces %3d/%3d shared %3d  gap %9.2e  %s\n", ...
            name, family, size(fm.elements, 2), size(sm.elements, 2), nF, nS, ...
            info.n_shared, info.area_gap, ternary(ok_shared && ok_gap, "ok", "BAD"));
        if ~(ok_shared && ok_gap), n_fail = n_fail + 1; end
        if family == "fitted", fm_fitted = fm; end
    end
    assert(~isempty(fm_fitted), "test_interface_meshes:noFitted", "the fitted family did not run");
    poly = interface_polygon(fm_fitted, curve);
    if poly.closed ~= (name ~= "curved")
        fprintf(2, "  interface_polygon closed = %d, expected %d\n", poly.closed, name ~= "curved");
        n_fail = n_fail + 1;
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_interface_meshes:failed", "%d checks failed", n_fail); end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
