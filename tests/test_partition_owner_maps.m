% Per-side point-to-interface maps (get_partitions mapping_1 / mapping_2).
% Each side's composite map x_i -> Sigma is affine in t on its faces, so
% its inverse on a patch is the affine s_i(t). Round trip at every patch
% node, from the patch endpoint data: the side-1 node x_1 -> x_bar =
% P_1(x_1) -> x_2 = P_2^{-1}(x_bar) -> P_2(x_2) -> P_1^{-1} must return
% x_1, and x_bar, x_2 must be the partitions' own interface and side-2
% points (the side-1 and side-2 rules pair the same points); P_1 must
% agree with map_mesh_to_interface's vertex rule up to the breakpoint
% merge (vertex_tol).
% Cases: gapped circle with vertices off Sigma, closest/normal both ways;
% gapped cosine graph with mapping_1 = "vertical", mapping_2 = "closest".
% With mapping_1 = mapping_2 = the curve's own rule the partitions equal
% the default call's bit for bit.
clc; clear;
n_fail = 0;
tol = 1e-12;
circle = CircularCurve([.5; .5], .4, mapping = "closest");
sine = SmoothCurveHorizontal("0.5 + .05*cos(2*pi*t)", mapping = "closest");
cases = {};
for h = [1/16, 1/32]
    [fm, sm] = interface_meshes(circle, "fitted", h, gap = 0.2*h, fluid_side = "outer");
    cases(end + 1, :) = {sprintf("circle h = 1/%d", round(1/h)), circle, fm, sm, ...
        ["closest", "normal"; "normal", "closest"]}; %#ok<SAGROW>
    % the graph is level at both walls, so the wall vertices map to the
    % ends under either rule (equal range)
    [fm, sm] = interface_meshes(sine, "fitted", h, gap = 0.2*h, fluid_side = "upper");
    cases(end + 1, :) = {sprintf("cosine h = 1/%d", round(1/h)), sine, fm, sm, ...
        ["vertical", "closest"]}; %#ok<SAGROW>
end
for c = 1:size(cases, 1)
    [name, curve, fm, sm, maps] = cases{c, :};
    for m = 1:size(maps, 1)
        [p1, p2] = curve.get_partitions(fm, sm, quadrature_patch = 3, ...
            mapping_1 = maps(m, 1), mapping_2 = maps(m, 2));
        % the patch composite maps: s_1 <-> t <-> s_2, affine per patch
        n_q = size(p1.quadrature.local, 2);
        e = repelem(1:size(p1.faces, 3), n_q);
        [v1, v2] = deal(reshape(p1.vertices, 2, []), reshape(p2.vertices, 2, []));
        tv = reshape(p1.vertices_interface, 2, []);
        affine = @(x, from, to) to(1, e) + (x - from(1, e))./(from(2, e) - from(1, e)).*(to(2, e) - to(1, e));
        X1 = reshape(p1.quadrature.physical, 2, []);
        s1 = reshape(p1.quadrature.local, 1, []);
        t = affine(s1, v1, tv);                       % x_bar = P_1(x_1)
        s2 = affine(t, tv, v2);                       % x_2 = P_2^{-1}(x_bar)
        X2 = face_point(sm, p2.faces(e), s2);
        X1_back = face_point(fm, p1.faces(e), affine(affine(s2, v2, tv), tv, v1));
        % the vertex rules themselves, before the breakpoint merge snaps
        % near-coincident vertex images together (vertex_tol)
        seg_1 = curve.map_mesh_to_interface(fm, mapping = maps(m, 1));
        err = [max(vecnorm(curve.get_2d_points(t) - reshape(p2.quadrature.interface, 2, []))), ...
               max(vecnorm(X2 - reshape(p2.quadrature.physical, 2, []))), ...
               max(vecnorm(X1_back - X1))];
        snap = max(abs(t - to_t(seg_1, p1.faces(e), s1)));
        fprintf("%-15s %-8s/%-8s  |x_bar - rule| %.1e  |x_2 - rule 2| %.1e  round trip %.1e  snap %.1e\n", ...
            name, maps(m, 1), maps(m, 2), err, snap);
        if snap > 1e-7, n_fail = n_fail + 1; end
        if any(err > tol), n_fail = n_fail + 1; end
    end
    [q1, q2] = curve.get_partitions(fm, sm);
    [r1, r2] = curve.get_partitions(fm, sm, mapping_1 = curve.mapping, mapping_2 = curve.mapping);
    if ~isequal(q1, r1) || ~isequal(q2, r2)
        fprintf(2, "  %s: explicit own maps differ from the default\n", name);
        n_fail = n_fail + 1;
    end
end
% a curve without a mapping property takes no override
try
    LevelSetCurve(@(x, y) x.^2 + y.^2 - 1).map_mesh_to_interface(cases{1, 3}, mapping = "closest");
    n_fail = n_fail + 1;
catch err
    if err.identifier ~= "AbstractCurve:noMapping", n_fail = n_fail + 1; end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_partition_owner_maps:failed", "%d checks failed", n_fail); end

function t = to_t(seg, f, s)
% composite map, face-local coordinate -> curve parameter (affine per segment)
t = zeros(size(s));
for i = 1:numel(s)
    q = find(seg.face == f(i));
    lo = min(seg.s(:, q), [], 1);  hi = max(seg.s(:, q), [], 1);
    q = q(find(s(i) >= lo - 1e-14 & s(i) <= hi + 1e-14, 1));
    t(i) = seg.t(1, q) + (s(i) - seg.s(1, q)) / diff(seg.s(:, q)) * diff(seg.t(:, q));
end
end

function X = face_point(mesh, f, s)
V = mesh.faces(1:2, f);
X = mesh.coordinates(1:2, V(1, :)).*(1 - s) + mesh.coordinates(1:2, V(2, :)).*s;
end
