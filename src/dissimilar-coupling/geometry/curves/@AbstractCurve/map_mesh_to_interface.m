function segments = map_mesh_to_interface(obj, mesh, options)
%-------------------------------------------------------------------------
%Description: Maps each interface face of `mesh` onto the arc of the
%             interface that it covers.
%
%             Returns parallel arrays, one column per segment:
%               .face  1 x N  mesh face index
%               .t     2 x N  [t_lo; t_hi], the arc covered, t_lo < t_hi
%               .s     2 x N  face-local coordinates matching t_lo and
%                             t_hi. s runs 0 -> 1 from the face's first
%                             vertex to its second, so s_lo > s_hi is
%                             legal: that is how direction is recorded,
%                             in place of a separate orientation flag.
%
%             A face contributes one segment, except on a periodic
%             interface where a face straddling the seam contributes two,
%             split there. Face-local coordinates are affine in the
%             interface parameter across a segment.
%
%             The vertex parameters come from map_points_to_interface;
%             under mapping = "normal" each vertex is cast along the
%             normalized average of its incident faces' unit normals.
%             The option `mapping` overrides the curve's own rule for this
%             mesh ("" keeps it; curves without a mapping property take
%             only "").
%
%             Side-constrained projection (periodic interfaces): the raw
%             per-vertex closest point can grab the wrong branch of the
%             curve where two stretches pass close together (a cusped tail
%             thinner than the mesh, a concave lobe), and the mapped
%             intervals then fold (get_partitions:foldedProjection). Along
%             each chain of interface faces the correct map is cyclically
%             monotone in t, so the repair keeps the longest monotone
%             subsequence of vertex parameters and re-projects the others
%             within the window between their kept neighbours. Faces
%             bridging a truncated stretch (the end cap of an unfitted
%             mesh at an under-resolved tail) are thereby stretched
%             monotonically across it, through the seam if needed.
%             Already-monotone mappings are untouched.
%
%             Ties (tip_mapping option, periodic interfaces). Where
%             several consecutive interface vertices share one closest
%             point (a sharp tip or cusp of Sigma, a cap wall
%             perpendicular to it), the faces between them have a
%             zero-length image:
%               "spread"   (default) the tied parameters are spread
%                          across the bracket to the next vertex, so
%                          every face keeps an arc of positive length
%                          (see spread_ties)
%               "collapse" the tied parameters are snapped to one value
%                          t*, also across the parameter seam, and each
%                          face between two tied vertices is returned as
%                          a collapsed segment: its whole length is sent
%                          to the single point gamma(t*)
%               .collapsed_face  1 x M  faces sent to a single point
%               .collapsed_t     1 x M  their parameter t*
%             (both empty under "spread"). "collapse" is for periodic
%             interfaces only and is refused on open ones: there a tie is
%             a run of faces beyond an end of Sigma, which no point of
%             Sigma pairs with (spread reports them as orphanFace).
%-------------------------------------------------------------------------
arguments
    obj AbstractCurve
    mesh SimplicialGrid
    options.tip_mapping (1,1) string ...
        {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
    options.mapping (1,1) string = ""
end
collapse = options.tip_mapping == "collapse";
if collapse && ~is_periodic(obj)
    error("map_mesh_to_interface:collapseOpenCurve", ...
        "tip_mapping = ""collapse"" needs a periodic (closed) interface; " + ...
        "%s is open.", class(obj));
end
if options.mapping ~= ""
    if ~isprop(obj, "mapping")
        error("AbstractCurve:noMapping", ...
            "%s has a single point-to-interface map; mapping = ""%s"" cannot apply.", ...
            class(obj), options.mapping);
    end
    obj.mapping = options.mapping;
end
% mapping = "normal": a vertex's ray must meet Sigma within this many of
% its longest incident face lengths
NORMAL_REACH = 2;
inter_faces = mesh.intfaces;
V = mesh.faces(1:2, inter_faces);
[vertex_ids, ~, local] = unique(V(:));
e = reshape(local, 2, []);                 % faces in chain-local vertex ids
X = mesh.coordinates(1:2, vertex_ids);
if isprop(obj, "mapping") && obj.mapping == "normal"
    [N, h] = vertex_normals(mesh, inter_faces, e, numel(vertex_ids));
    v_t = reshape(map_points_to_interface(obj, X, normals = N, ...
        reach = NORMAL_REACH * h), 1, []);
else
    v_t = reshape(map_points_to_interface(obj, X), 1, []);
end

if is_periodic(obj)
    v_t = unfold_vertices(obj, e, X, v_t, collapse);
end

t_start = v_t(e(1, :));   % s = 0
t_end   = v_t(e(2, :));   % s = 1

t0 = obj.t_range(1);
period = diff(obj.t_range);
if is_periodic(obj)
    % Take the shorter of the two arcs joining the endpoints,
    % unwrapping past the seam where that is the short way
    % round. Afterwards |t_end - t_start| <= period/2, so a
    % face can straddle at most one seam.
    delta = t_end - t_start;
    wraps = abs(delta) > period/2;
    delta(wraps) = delta(wraps) - sign(delta(wraps)) * period;
    t_end = t_start + delta;
end

n_faces = numel(inter_faces);
face = zeros(1, 2*n_faces);
[tt, ss] = deal(zeros(2, 2*n_faces));
n = 0;
[collapsed_face, collapsed_t] = deal(zeros(1, 0));
for i = 1:n_faces
    a = t_start(i);
    b = t_end(i);
    if a == b   % zero-length projection
        if collapse
            collapsed_face(end + 1) = inter_faces(i); %#ok<AGROW>
            collapsed_t(end + 1) = t0 + mod(a - t0, period); %#ok<AGROW>
        end
        continue
    end

    % Seam splitting is a periodic-only construction: on an
    % open interface find_closest_point can return endpoints a
    % few eps outside t_range, and unwrapping those would emit
    % a breakpoint a whole period away.
    if is_periodic(obj) && b < t0   % wrapped below the seam
        s_seam = (t0 - a) / (b - a);
        pieces = [a, t0, 0, s_seam; t0 + period, b + period, s_seam, 1];
    elseif is_periodic(obj) && b > t0 + period  % wrapped above
        s_seam = (t0 + period - a) / (b - a);
        pieces = [a, t0 + period, 0, s_seam; t0, b - period, s_seam, 1];
    else
        pieces = [a, b, 0, 1];
    end

    for q = 1:size(pieces, 1)
        ta = pieces(q, 1); tb = pieces(q, 2);
        sa = pieces(q, 3); sb = pieces(q, 4);
        if ta > tb
            [ta, tb, sa, sb] = deal(tb, ta, sb, sa);
        end
        if ta == tb, continue; end
        n = n + 1;
        face(n) = inter_faces(i);
        tt(:, n) = [ta; tb];
        ss(:, n) = [sa; sb];
    end
end

segments.face = face(1:n);
segments.t = tt(:, 1:n);
segments.s = ss(:, 1:n);
segments.collapsed_face = collapsed_face;
segments.collapsed_t = collapsed_t;
end

%=========================================================================%
function [N, h] = vertex_normals(mesh, inter_faces, e, n_vertices)
% Per interface vertex: the normalized average of the outward unit normals
% of its incident interface faces (outward from `mesh`, so consistently
% oriented along every chain), and the longest incident face length.
[K, side] = mesh.get_element_by_face(inter_faces);
n_face = reshape(mesh.get_normal_by_element(K, side).', 2, []);   % scaled by length
len = vecnorm(n_face);
n_face = n_face ./ len;
N = [accumarray(e(:), reshape(n_face([1 1], :), [], 1), [n_vertices 1]), ...
     accumarray(e(:), reshape(n_face([2 2], :), [], 1), [n_vertices 1])].';
N = N ./ vecnorm(N);
h = accumarray(e(:), reshape(len([1 1], :), [], 1), [n_vertices 1], @max).';
end

function v_t = unfold_vertices(obj, e, X, v_t, collapse)
% Cyclic-monotone repair of the per-vertex parameters (see the doc
% comment above). collapse: ties are snapped (snap_ties) instead of
% spread (spread_ties). Works in relative coordinates from a cut vertex,
% v = mod(t - t_cut, period): a once-winding monotone map is exactly a
% non-decreasing v along the chain, with no unwrap step whose winding
% folded data could destroy. The cut is chosen to maximize the kept
% (already monotone) subsequence. No-op when a chain is already
% monotone from some cut.
t0 = obj.t_range(1);
period = diff(obj.t_range);
opts = optimset('TolX', 1e-12, 'MaxIter', 2e3, 'MaxFunEvals', 2e3);

[chains, is_cycle] = trace_chains(e, numel(v_t));
for c = 1:numel(chains)
    chain = chains{c};
    n = numel(chain);
    if n < 4, continue; end
    t = v_t(chain);

    % Candidate cuts: for a cycle any rotation, in either orientation,
    % since the walk direction of trace_chains is arbitrary; for an open
    % chain only the two orientations. The seam-adjacent vertex is tried
    % first, where an already-monotone map is confirmed in one try.
    best_len = -1;
    [~, seam] = min(mod(t - t0, period));
    for direction = 1:2
        if is_cycle(c)
            cuts = [seam, 1:n];
        else
            cuts = 1;
        end
        for cut = cuts
            if direction == 1
                order = [cut:n, 1:cut - 1];
            elseif is_cycle(c)
                order = [cut:-1:1, n:-1:cut + 1];
            else
                order = n:-1:1;
            end
            v = mod(t(order) - t(order(1)), period);
            len = monotone_length(v);
            if len > best_len
                best_len = len; best_order = order; best_v = v;
            end
            if len == n, break; end
        end
        if best_len == n, break; end
    end
    order = best_order;
    v = best_v;
    base = t(order(1));

    if best_len == n
        % Monotone from some cut: untouched, except that (near-)zero-length
        % runs from perpendicular cap faces are spread (see spread_ties).
        % Only moved vertices are written back, so a clean chain stays
        % bit-for-bit identical.
        if collapse
            % copy the representative's parameter bit for bit
            [~, moved, source] = snap_ties(v, period, is_cycle(c));
            v_t(chain(order(moved))) = v_t(chain(order(source(moved))));
            continue
        end
        [v, moved] = spread_ties(v, period);
        if any(moved)
            idx = chain(order(moved));
            v_t(idx) = t0 + mod(base + v(moved) - t0, period);
        end
        continue
    end

    keep = monotone_mask(v);
    kept = find(keep);
    for j = find(~keep)
        % Bracket by the nearest kept neighbours in relative
        % coordinates; the chain start and end close the window at 0
        % and period.
        before = kept(find(kept < j, 1, "last"));
        after  = kept(find(kept > j, 1, "first"));
        if isempty(before); lo = 0;      else; lo = v(before); end
        if isempty(after);  hi = period; else; hi = v(after);  end
        if hi - lo > 0
            x = X(:, chain(order(j)));
            distance = @(w) sum((x - obj.get_2d_points(base + w)).^2);
            v(j) = fminbnd(distance, lo, hi, opts);
        else
            v(j) = lo;   % degenerate window: collapse (face drops later)
        end
        % Sequential monotonicity inside a run of repaired vertices
        keep(j) = true;
        kept = find(keep);
    end
    if collapse
        [v, moved, source] = snap_ties(v, period, is_cycle(c));
        v_t(chain(order)) = t0 + mod(base + v - t0, period);
        v_t(chain(order(moved))) = v_t(chain(order(source(moved))));
    else
        v = spread_ties(v, period);
        v_t(chain(order)) = t0 + mod(base + v - t0, period);
    end
end
end

function [v, moved] = spread_ties(v, period)
% Faces perpendicular to the interface (cap walls bridging a truncated
% stretch) project to runs of (nearly) equal parameters; their faces
% would drop or collapse under get_partitions' breakpoint merge and leave
% their trace dofs unconstrained (orphanFace). Spread each run across its
% bracket; tol sits well above the merge tolerance and well below any
% resolved face's arc. `moved` marks the changed entries (a run's first
% entry stays put).
n = numel(v);
moved = false(1, n);
tol = 1e-6 * period;
j = 1;
while j <= n
    k = j;
    while k < n && v(k + 1) - v(j) < tol
        k = k + 1;
    end
    if k > j
        if k < n; run_hi = v(k + 1); else; run_hi = period; end
        v(j:k) = v(j) + (0:(k - j)) * (run_hi - v(j)) / (k - j + 1);
        moved(j + 1:k) = true;
    end
    j = k + 1;
end
end

function [v, moved, source] = snap_ties(v, period, is_cycle)
% The collapse alternative to spread_ties: each run of (nearly) equal
% parameters, with the same tol, is snapped to its first entry, so the
% faces inside the run have a zero-length image and map_mesh_to_interface
% returns them as collapsed segments. On a cycle (is_cycle) the run may
% straddle the cut: entries within tol below the period tie with
% v(1) = 0, the same point of Sigma, and are snapped to it; an open chain
% does not wrap. source(i) is the entry whose
% value entry i takes (the caller copies parameters by it, so tied
% vertices share their parameter bit for bit); moved marks the entries
% with source(i) ~= i.
n = numel(v);
source = 1:n;
tol = 1e-6 * period;
j = 1;
while j <= n
    k = j;
    while k < n && v(k + 1) - v(j) < tol
        k = k + 1;
    end
    source(j + 1:k) = j;
    v(j + 1:k) = v(j);
    j = k + 1;
end
wrap = find(period - v < tol & source ~= 1);
if is_cycle && ~isempty(wrap) && v(1) == 0
    source(wrap) = 1;
    v(wrap) = 0;
end
moved = source ~= (1:n);
end

function len = monotone_length(v)
% Length of the longest non-decreasing subsequence (patience tails only)
n = numel(v);
tails = zeros(1, n);
len = 0;
for i = 1:n
    lo = 1; hi = len;
    while lo <= hi
        mid = floor((lo + hi) / 2);
        if tails(mid) <= v(i); lo = mid + 1; else; hi = mid - 1; end
    end
    tails(lo) = v(i);
    len = max(len, lo);
end
end

function keep = monotone_mask(v)
% Longest non-decreasing subsequence (patience sorting, O(n log n)):
% the maximal set of vertices already consistent with a monotone map.
n = numel(v);
[tails, tails_at, prev] = deal(zeros(1, n));
len = 0;
for i = 1:n
    lo = 1; hi = len;
    while lo <= hi
        mid = floor((lo + hi) / 2);
        if tails(mid) <= v(i); lo = mid + 1; else; hi = mid - 1; end
    end
    tails(lo) = v(i);
    tails_at(lo) = i;
    if lo > 1; prev(i) = tails_at(lo - 1); end
    len = max(len, lo);
end
keep = false(1, n);
i = tails_at(len);
while i > 0
    keep(i) = true;
    i = prev(i);
end
end

function [chains, is_cycle] = trace_chains(e, n_vertices)
% Orders the interface vertices along their face chains. Each vertex of
% a manifold chain touches at most two faces; a vertex touching more (a
% pinched staircase) breaks the walk there and bounds the chains.
n_faces = size(e, 2);
degree = zeros(n_vertices, 1);
incident = zeros(n_vertices, 2);
for f = 1:n_faces
    for side = 1:2
        v = e(side, f);
        degree(v) = degree(v) + 1;
        if degree(v) <= 2
            incident(v, degree(v)) = f;
        end
    end
end

visited = false(1, n_faces);
chains = {};
is_cycle = logical.empty(1, 0);
% Open chains first (endpoints of degree ~= 2), then leftover cycles
starts = [reshape(find(degree(e(1, :)) ~= 2 | degree(e(2, :)) ~= 2), 1, []), 1:n_faces];
for f0 = starts
    if visited(f0), continue; end
    v0 = e(1, f0);
    if degree(e(2, f0)) ~= 2 && degree(v0) == 2
        v0 = e(2, f0);   % start the walk at the chain's end
    end
    seq = v0;
    f = f0;
    closed = false;
    while true
        visited(f) = true;
        v = e(1, f) + e(2, f) - seq(end);   % the face's other vertex
        if v == seq(1)
            closed = true;
            break
        end
        seq(end + 1) = v; %#ok<AGROW>
        if degree(v) ~= 2, break; end
        f_next = incident(v, 1) + incident(v, 2) - f;
        if f_next == 0 || visited(f_next), break; end
        f = f_next;
    end
    chains{end + 1} = seq; %#ok<AGROW>
    is_cycle(end + 1) = closed; %#ok<AGROW>
end
end
