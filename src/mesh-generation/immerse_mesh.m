function new_mesh = immerse_mesh(obj, phi)
%IMMERSE_MESH  Carve the domain {phi < 0} out of a background grid.
%
%   new_mesh = immerse_mesh(obj, phi) keeps the elements of the background
%   SimplicialGrid obj whose three vertices all lie strictly inside
%   {phi < 0}; a vertex with phi == 0 does not count as inside. phi is a
%   function handle of one argument, called as phi(X) on the 2 x 3 x n
%   array of the vertex coordinates of the n elements (first dimension: x
%   and y) and returning the values with first dimension 1, or a
%   LevelSetCurve (evaluated through its evaluate method).
%
%   new_mesh is the SimplicialGrid of the kept elements, as extracted by
%   obj(elements). Its faces on the boundary of the kept region, whatever
%   their previous tag, are tagged "Dirichlet"; interior faces stay
%   untagged.
%
%See also: LevelSetCurve, TransferEssentialBoundaryConditions, SimplicialGrid/split
    arguments
        obj SimplicialGrid
        phi {mustBeA(phi, ["function_handle", "LevelSetCurve"])}
    end
    if isa(phi, "LevelSetCurve")
        phi = @(X) phi.evaluate(X);
    end

    vertices = mat2cell(obj.vertices, [1; 1]);
    vertices = cat(1, vertices{:});
    z_vertices = phi(vertices);
    in_vertices = (z_vertices <= 0);
    on_vertices = (z_vertices == 0);
    in_vertices = in_vertices & ~on_vertices;

    int_K = find(all(in_vertices));

    new_mesh = obj(int_K);
    new_mesh.faces(3, :) = sign(new_mesh.faces(3, :));
end