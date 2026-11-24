function export_fields_csv(stem, fluid_solution, solid_solution, ...
    fluid_mesh, solid_mesh, options)
%-------------------------------------------------------------------------
%Description: Writes the coupled solution's fields as plain CSV, one row
%             per element:
%
%               [x1 y1 x2 y2 x3 y3, f1 f2 f3, g1 g2 g3, ux1..3, uy1..3]
%
%             the element's corners, two scalar fields at those corners --
%             fluid: f = |u_F|, g = p_F; solid: f = |u_S|, g = von Mises --
%             and then the vector components, for streamlines and quivers.
%             Two files are written, <stem>_fluid.csv and <stem>_solid.csv.
%             Intended for headless runs, where MATLAB cannot create
%             figures.
%
%             The fields are evaluated, never read out of raw dof rows:
%             the PEERS stress array is not four equal blocks (it carries
%             the curl-curl bubble for k >= 1).
%
%See also: extend_fields_along_paths (the TPM extension across the gap)
%-------------------------------------------------------------------------
arguments
    stem            (1,1)   string          % path prefix; _fluid/_solid.csv appended
    fluid_solution  (1,1)   FESolution
    solid_solution  (1,1)   FESolution
    fluid_mesh              SimplicialGrid
    solid_mesh              SimplicialGrid
    options.poisson_ratio (1,1) double = 0.3   % for the plane-strain sigma_33
end
at = @(F, i, j) reshape(F(:, :, i, j, :), 3, []);            % 3 x n_elements
vertex = @(field) field.evaluate(reference_element(2), coordinate_system = "vertex");
corners = @(mesh) reshape(mesh.coordinates(:, mesh.elements(1:3, :)), 6, []);

U = vertex(fluid_solution.Unknowns.field(2));
S = vertex(fluid_solution.Unknowns.field(1));
ux = at(U, 1, 1);  uy = at(U, 2, 1);
speed = hypot(ux, uy);
pressure = -0.5*(at(S, 1, 1) + at(S, 2, 2));
writematrix([corners(fluid_mesh); speed; pressure; ux; uy].', stem + "_fluid.csv");

U = vertex(solid_solution.Unknowns.field(2));
S = vertex(solid_solution.Unknowns.field(1));
s11 = at(S, 1, 1);  s22 = at(S, 2, 2);  s12 = 0.5*(at(S, 1, 2) + at(S, 2, 1));
s33 = options.poisson_ratio*(s11 + s22);
von_mises = sqrt(0.5*((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2) + 3*s12.^2);
ux = at(U, 1, 1);  uy = at(U, 2, 1);
writematrix([corners(solid_mesh); hypot(ux, uy); von_mises; ux; uy].', ...
    stem + "_solid.csv");
end
