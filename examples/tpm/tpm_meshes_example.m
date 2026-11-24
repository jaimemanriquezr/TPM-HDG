% Mesh generation gallery: the unit square, an L-shaped polygon with two
% boundary markers, and a square with a square hole. Plots each mesh.
clc; clear;
mesh_software = "mesh2d";
h = .5/2;

% Square
mesh_square = UnitSquareMesh(h);
figure();
plot(mesh_square);

% Polygon with two boundaries
boundary_markers = [1, 1, 2, 2, 1, 1];
mesh_L = LShapeGrid(h, boundary_markers);
figure();
plot(mesh_L);

% Domain with a hole (mesh2d only)
big_square = [-2 2 2 -2; -2 -2 2 2];
small_square = [-1 1 1 -1; -1 -1 1 1];
connectivity  = [[1; 2], [2; 3], [3; 4], [4; 1], [5; 6], [6; 7], [7; 8], [8; 5]];
boundary_markers = [1, 1, 1, 1, 2, 2, 2, 2];
mesh_hole = SimplicialGrid([big_square, small_square], [connectivity; boundary_markers], ...
    element_diameter = h, ...
    mesh_generator = mesh_software);
figure();
plot(mesh_hole);
