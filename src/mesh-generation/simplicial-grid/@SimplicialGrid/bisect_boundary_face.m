function obj = bisect_boundary_face(obj, faces, at)
%-------------------------------------------------------------------------
%Description: Splits each boundary face in `faces` at a point p, and its
%             element into two along p -> the opposite vertex. The mesh is
%             rebuilt ("array"): vertex numbers survive (p is appended),
%             face numbers do not. Every surviving face keeps its tag, the
%             two halves take the split face's tag, the new interior face
%             is untagged. tau is reset, so call this before
%             add_stabilization.
%
%             The split point is given per face either as the face-local
%             coordinate s in (0, 1) (s = 0 at faces(1, f)), or as a
%             function g(X) (2 x n -> 1 x n) whose sign change along the
%             face is p. Two common choices:
%               equidistant  @(X) Ci.distance(X) - Cj.distance(X)
%               parameter    @(X) C.map_points_to_interface(X) - t_target
%
%Usage:       mesh = bisect_boundary_face(mesh, f, 0.5)
%             mesh = bisect_boundary_face(mesh, [f1 f2], {g, g})
%
%See also: tag_by_curves, uniform_refine
%-------------------------------------------------------------------------
arguments
    obj SimplicialGrid
    faces (1, :) double
    at                      % (1, :) double, or (1, :) cell of function handles
end
assert(numel(at) == numel(faces), "SimplicialGrid:bisect_boundary_face:size", ...
    "one split point per face");
% the faces are identified by their vertex pairs, which survive the rebuilds
pairs = obj.faces(1:2, faces);
for i = 1:numel(faces)
    obj = bisect_one(obj, pairs(:, i), at(i));
end
end

%=========================================================================%
function obj = bisect_one(obj, pair, at)
n_v = size(obj.coordinates, 2);
key = @(F) min(F(1:2, :), [], 1) + (n_v + 1)*max(F(1:2, :), [], 1);
face = find(key(obj.faces) == key(pair));
assert(isscalar(face), "SimplicialGrid:bisect_boundary_face:missing", ...
    "face (%d, %d) is not in the mesh", pair(1), pair(2));
[K, ~] = get_element_by_face(obj, face);
assert(isscalar(K) && ~iscell(K), "SimplicialGrid:bisect_boundary_face:interior", ...
    "face (%d, %d) is not a boundary face", pair(1), pair(2));

ab = obj.faces(1:2, face);
a = obj.coordinates(:, ab(1));  b = obj.coordinates(:, ab(2));
if iscell(at)
    g = at{1};
    ga = g(a);  gb = g(b);
    assert(sign(ga)*sign(gb) < 0, "SimplicialGrid:bisect_boundary_face:noCrossing", ...
        "the split function does not change sign along face (%d, %d)", ab(1), ab(2));
    s = fzero(@(s) g(a + s*(b - a)), [0, 1], optimset("TolX", 1e-15));
else
    s = at;
end
assert(s > 1e-6 && s < 1 - 1e-6, "SimplicialGrid:bisect_boundary_face:noCrossing", ...
    "split point s = %g is not inside face (%d, %d)", s, ab(1), ab(2));
p = a + s*(b - a);
n_p = n_v + 1;
coordinates = [obj.coordinates, p];

% p takes one face vertex's slot in each half: orientation is kept
elements = obj.elements;
E = elements(:, K);
E1 = E;  E1(E == ab(2)) = n_p;
E2 = E;  E2(E == ab(1)) = n_p;
elements(:, K) = E1;
elements(:, end + 1) = E2;

new = SimplicialGrid(coordinates, elements, mesh_generator = "array");
new.boundary_tags = obj.boundary_tags;
% tags by vertex pair; the halves inherit the split face's tag
key = @(F) min(F(1:2, :), [], 1) + (n_p + 1)*max(F(1:2, :), [], 1);
tag = obj.faces(3, face);
old_faces = [obj.faces, [ab(1); n_p; tag], [n_p; ab(2); tag]];
[found, loc] = ismember(key(new.faces), key(old_faces));
assert(nnz(~found) == 1, "SimplicialGrid:bisect_boundary_face:faces", ...
    "expected exactly one new interior face");
new.faces(3, :) = 0;
new.faces(3, found) = old_faces(3, loc(found));
obj = new;
end
