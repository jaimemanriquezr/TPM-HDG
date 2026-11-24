function obj = discretize(obj, options)
% Induced parametrization of the zero contour: an ordered polyline with
% an arc-length coordinate s in [0, L], serving the ParametrizedCurve
% query surface (t_range, get_2d_points, tangent,
% map_points_to_interface). Vertices are Newton-projected onto
% {phi = 0}, so they lie on the curve, not on the contour-grid chords.
% Scope: a single closed component inside obj.box.
arguments
    obj LevelSetCurve
    options.n_grid (1,1) double = 512   % contour-extraction grid resolution
    % optional seam anchor: roll the parametrization so s = 0 is the
    % curve point nearest this point (contourc otherwise starts the
    % polyline at an arbitrary place, which moves the periodic seam)
    options.anchor double = double.empty(2, 0)
end
b = obj.box;
xv = linspace(b(1), b(2), options.n_grid);
yv = linspace(b(3), b(4), options.n_grid);
[Xg, Yg] = meshgrid(xv, yv);
Z = reshape(obj.evaluate([Xg(:).'; Yg(:).']), size(Xg));
C = contourc(xv, yv, Z, [0 0]);

% contourc column format -> component list
comps = {};
i = 1;
while i < size(C, 2)
    n = C(2, i);
    comps{end+1} = C(:, i+1 : i+n); %#ok<AGROW>
    i = i + n + 1;
end
assert(~isempty(comps), "LevelSetCurve:discretize", ...
    "the zero contour does not intersect box %s", mat2str(b));
assert(isscalar(comps), "LevelSetCurve:discretize", ...
    ["the zero contour has %d components in box %s; only a single " ...
     "closed component is supported"], numel(comps), mat2str(b));
P = comps{1};
chord = @(Q) sum(vecnorm(diff(Q, 1, 2)));
assert(norm(P(:, 1) - P(:, end)) < 1e-8 * chord(P), ...
    "LevelSetCurve:discretize", ...
    "open zero contours are not supported (component endpoints differ)");
P = P(:, 1:end-1);            % periodic storage: closing vertex implicit

P = obj.project(P);           % vertices exactly on {phi = 0}
% projection can collapse neighbors; drop near-duplicates
gaps = vecnorm(diff([P, P(:, 1)], 1, 2));
P = P(:, gaps > 1e-10 * max(gaps));

% canonical orientation, matching the parametrized-lineage convention:
% the outward (grad phi) normal is the tangent rotated by -90 degrees
T1 = P(:, 2) - P(:, 1);
N1 = obj.normal_at(P(:, 1));
if [T1(2); -T1(1)].' * N1 < 0
    P = [P(:, 1), fliplr(P(:, 2:end))];
end

if ~isempty(options.anchor)
    [~, j] = min(vecnorm(P - obj.project(options.anchor(:))));
    P = circshift(P, -(j-1), 2);
end

seg = vecnorm(diff([P, P(:, 1)], 1, 2));
obj.polyline = P;
obj.arc = [0, cumsum(seg)];   % (1, m+1); arc(end) = L
end
