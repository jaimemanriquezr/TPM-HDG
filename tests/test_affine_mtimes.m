%-------------------------------------------------------------------------
%Description: SimplicialGrid mtimes against a rebuild from arrays.
%
%             A * T must be the mesh one gets by constructing a grid from
%             the mapped vertices, SimplicialGrid(pagemtimes(A, T.coordinates),
%             T.elements, mesh_generator = "array"): same connectivity and
%             the same volume, face lengths and (scaled) normals. The face
%             markers are not compared: the array constructor does not know
%             T's tags, mtimes keeps them.
%-------------------------------------------------------------------------
n_failed = 0;

T = StructuredSimplicialGrid(3, 2);

% orientation-preserving maps only: the volume is always positive
maps = {[1 .3; 0 .5], [1 0; 0 .5], [2 0; 0 2]};
names = ["skew", "fitted_split scaling", "isotropic scaling"];

for i = 1:numel(maps)
    A = maps{i};
    fprintf("A = %s\n", names(i));
    lhs = A * T;
    ref = SimplicialGrid(pagemtimes(A, T.coordinates), T.elements, mesh_generator = "array");

    n_failed = n_failed + check(names(i) + ": coordinates", ...
        max(abs(lhs.coordinates - ref.coordinates), [], "all") < 1e-14);
    n_failed = n_failed + check(names(i) + ": elements", isequal(lhs.elements, ref.elements));
    n_failed = n_failed + check(names(i) + ": faces (vertices)", ...
        isequal(lhs.faces(1:2, :), ref.faces(1:2, :)));
    n_failed = n_failed + check(names(i) + ": boundaries", isequal(lhs.boundaries, ref.boundaries));
    n_failed = n_failed + check(names(i) + ": incidence", isequal(lhs.incidence, ref.incidence));
    n_failed = n_failed + check(names(i) + ": orientation", isequal(lhs.orientation, ref.orientation));
    n_failed = n_failed + check(names(i) + ": permutation", isequal(lhs.permutation, ref.permutation));
    n_failed = n_failed + check(names(i) + ": volume", ...
        max(abs(lhs.volume - ref.volume)) < 1e-14);
    n_failed = n_failed + check(names(i) + ": area", ...
        max(abs(lhs.area - ref.area)) < 1e-14);
    n_failed = n_failed + check(names(i) + ": normals", ...
        max(abs(lhs.normals - ref.normals), [], "all") < 1e-14);
    n_failed = n_failed + check(names(i) + ": tau", isequal(lhs.tau, ref.tau));
end

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
