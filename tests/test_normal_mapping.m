% mapping = "normal": map_mesh_to_interface casts each interface vertex
% along the averaged unit normal of its incident faces. On gapped fitted
% meshes of a circle and of a sine graph, per mesh: every vertex lands
% within reach (else ParametrizedCurve:noNormalIntersection), the vertex
% parameters do not fold (on the circle the full get_partitions pair also
% builds), the transfer path x -> gamma(t) is parallel to the averaged
% normal, and its angle to the exact normal at gamma(t) decreases with h
% (on the circle the gapped polygons are regular and the angle vanishes).
% Gapped open pairs fail get_partitions:unequalRange at the box walls
% under "closest" as well, so the graph is checked per mesh only.
clc; clear;
n_fail = 0;
cases = {"circle", @(m) CircularCurve([.5; .5], .3, mapping = m), "outer";
         "sine",   @(m) SmoothCurveHorizontal("0.5 + .05*sin(2*pi*t)", mapping = m), "lower"};
hs = [1/8, 1/16, 1/32];
for c = 1:size(cases, 1)
    [name, make_curve, side] = cases{c, :};
    curve = make_curve("normal");
    angle = zeros(1, numel(hs));
    for l = 1:numel(hs)
        h = hs(l);
        [fm, sm] = interface_meshes(curve, "fitted", h, gap = 0.2*h, fluid_side = side);
        for mesh = {fm, sm}
            [t, X, N] = vertex_map(curve, mesh{1});
            % no fold: the segments' parameter intervals do not overlap
            seg = curve.map_mesh_to_interface(mesh{1});
            [lo, order] = sort(seg.t(1, :));
            hi = seg.t(2, order);
            if any(lo(2:end) < hi(1:end-1) - 1e-12)
                fprintf(2, "  %s h = 1/%d: folded segments\n", name, round(1/h));
                n_fail = n_fail + 1;
            end
            % open curve: a wall vertex whose ray lands past an end is
            % clamped there, off its ray
            inside = curve.is_periodic() | (t > curve.t_range(1) & t < curve.t_range(2));
            [t, X, N] = deal(t(inside), X(:, inside), N(:, inside));
            D = curve.get_2d_points(t) - X;
            r = vecnorm(D);
            off_ray = max(abs(D(1, :).*N(2, :) - D(2, :).*N(1, :)));
            if off_ray > 1e-10 || max(r) > 2*h
                fprintf(2, "  %s h = 1/%d: off the ray %g, max |r| %g\n", ...
                    name, round(1/h), off_ray, max(r));
                n_fail = n_fail + 1;
            end
            nu = curve.unit_normal(t);
            angle(l) = max(angle(l), max(asin(abs(N(1, :).*nu(2, :) - N(2, :).*nu(1, :)))));
        end
        if name == "circle"
            try
                curve.get_partitions(fm, sm);
            catch err
                fprintf(2, "  %s h = 1/%d: get_partitions: %s\n", name, round(1/h), err.message);
                n_fail = n_fail + 1;
            end
        end
    end
    rate = log2(angle(1:end-1)./angle(2:end));
    fprintf("%-7s max path angle to the exact normal: %s, rates %s\n", ...
        name, mat2str(angle, 3), mat2str(rate, 3));
    % a regular polygon's averaged normals are exact (the circle)
    if any(rate < 0.9) && max(angle) > 1e-10
        n_fail = n_fail + 1;
    end
end

% no normals: the rule cannot run and says so
curve = CircularCurve([.5; .5], .3, mapping = "normal");
try
    curve.map_points_to_interface([.5; .9]);
    n_fail = n_fail + 1;
catch err
    if err.identifier ~= "SmoothCurve:normalsRequired"
        n_fail = n_fail + 1;
    end
end
% a ray that misses within reach fails loudly
try
    curve.map_points_to_interface([.5; .5], normals = [1; 0], reach = 0.1);
    n_fail = n_fail + 1;
catch err
    if err.identifier ~= "ParametrizedCurve:noNormalIntersection"
        n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_normal_mapping:failed", "%d checks failed", n_fail); end

function [t, X, N] = vertex_map(curve, mesh)
% the vertex parameters as map_mesh_to_interface assigns them (segment
% ends at s = 0 or 1), with the vertices and their averaged normals
seg = curve.map_mesh_to_interface(mesh);
F = mesh.faces(1:2, seg.face);
[vid, t] = deal([]);
for side = 1:2
    s_end = side - 1;
    at = seg.s == s_end;
    for row = 1:2
        vid = [vid, F(side, at(row, :))]; %#ok<AGROW>
        t = [t, seg.t(row, at(row, :))];  %#ok<AGROW>
    end
end
[vid, first] = unique(vid);
t = t(first);
X = mesh.coordinates(1:2, vid);
% averaged outward unit normals of the incident interface faces
faces = unique(seg.face);
[K, sn] = mesh.get_element_by_face(faces);
n_face = reshape(mesh.get_normal_by_element(K, sn).', 2, []);
n_face = n_face ./ vecnorm(n_face);
[~, j] = ismember(mesh.faces(1:2, faces), vid);
N = [accumarray(j(:), reshape(n_face([1 1], :), [], 1), [numel(vid) 1]), ...
     accumarray(j(:), reshape(n_face([2 2], :), [], 1), [numel(vid) 1])].';
N = N ./ vecnorm(N);
end
