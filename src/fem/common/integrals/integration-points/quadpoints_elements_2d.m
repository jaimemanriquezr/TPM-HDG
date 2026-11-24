function [xy, rs, w, F] = quadpoints_elements_2d(mesh, order, element_indices, options)
%QUADPOINTS_ELEMENTS_2D  Quadrature points on the elements of a mesh.
%
%   [xy, rs, w, F] = quadpoints_elements_2d(mesh, order, element_indices, ...
%                                           output = "cell")
%
%   Maps the reference rule get_2d_reference_quadrature(order, mesh.ElementType) to the
%   chosen elements (default all) with their affine maps F. With output
%   "cell" xy and rs are split into {x; y} and {r; s} cells and F into a
%   2-by-2 cell; "array" returns them numeric.
arguments
    mesh AbstractGrid
    order double
    element_indices (:,1) double = 1:size(mesh.elements, 2);
    options.output = "cell"
end
[rs, w] = get_2d_reference_quadrature(order, mesh.ElementType);
[F, xy_0] = get_affine_map(mesh, output="array");
F = F(:, :, element_indices);
xy_0 = xy_0(:, element_indices);
xy = pagemtimes(F, rs) + reshape(xy_0, 2, 1, []);

if options.output == "cell"
    xy = {reshape(xy(1, :), [], length(element_indices));
        reshape(xy(2, :), [], length(element_indices))};
    rs = {rs(1, :); rs(2, :)};
    F = {F(1, 1, :), F(1, 2, :); F(2, 1, :), F(2, 2, :)};
end
end
