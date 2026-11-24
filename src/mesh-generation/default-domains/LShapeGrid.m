function [mesh, domain] = LShapeGrid(h, boundary_markers)
%LSHAPEGRID  Mesh of the L-shaped domain [-1, 1]^2 without [0, 1] x [0, 1].
%
%   [mesh, domain] = LShapeGrid(h, boundary_markers)
%
%   h                element diameter (default 0.1).
%   boundary_markers 1-by-6 tags of the six edges, in the order of the
%                    vertices (-1,-1), (1,-1), (1,0), (0,0), (0,1), (-1,1)
%                    (default all 1).
%   mesh             SimplicialGrid generated with mesh2d.
%   domain           {vertices, boundary_markers} that define it.
arguments
    h double = .1
    boundary_markers (1, 6) = [1, 1, 1, 1, 1, 1]
end
vertices = [[-1; -1], [1; -1], [1; 0], [0; 0], [0; 1], [-1; 1]];
domain = {vertices, boundary_markers};
mesh = SimplicialGrid(domain{:}, ...
    element_diameter = h, ...
    mesh_generator = "mesh2d");
end
