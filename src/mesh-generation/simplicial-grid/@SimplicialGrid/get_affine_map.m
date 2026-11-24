function [F, x3] = get_affine_map(obj, K, options)
%GET_AFFINE_MAP  Affine maps from the reference element to mesh elements.
%
%   [F, x3] = get_affine_map(mesh, K, inverse = false, output = "array")
%
%   For element K (default all elements) with vertices x_1, x_2, x_3, the
%   map from reference coordinates (r, s) is x = B*[r; s] + x_3 with
%   B = [x_1 - x_3, x_2 - x_3]. F is B (or inv(B) when inverse is true),
%   returned as a 2-by-2-by-numel(K) array for output = "array" or as a
%   2-by-2 cell of column vectors for output = "cell"; x3 is the 2-by-numel(K)
%   array of the points x_3, the image of the reference origin.
arguments
    obj SimplicialGrid
    K (:,1) = 1:size(obj.elements,2)
    options.inverse = false
    options.output = "array"
end

if isempty(obj)
    F = eye(2);
    x3 = [0.; 0.];
    return
end

x_1 = obj.coordinates(1, obj.elements(1, K));
x_2 = obj.coordinates(1, obj.elements(2, K));
x_3 = obj.coordinates(1, obj.elements(3, K));
B_11 = x_1 - x_3;
B_12 = x_2 - x_3;

y_1 = obj.coordinates(2, obj.elements(1, K));
y_2 = obj.coordinates(2, obj.elements(2, K));
y_3 = obj.coordinates(2, obj.elements(3, K));
B_21 = y_1 - y_3;
B_22 = y_2 - y_3;

x3 = [x_3; y_3];
if options.inverse
    detB = (B_11.*B_22 - B_12.*B_21);
    G_11 = B_22 ./ detB;
    G_21 = -(B_21 ./ detB);
    G_12 = -(B_12 ./ detB);
    G_22 = B_11 ./ detB;

    F = [G_11(:), G_21(:), G_12(:), G_22(:)];
else
    F = [B_11(:), B_21(:), B_12(:), B_22(:)];
end

switch options.output
    case "array"
        F = reshape(F.', 2, 2, []);
    case "cell"
        F = reshape(mat2cell(F(:), length(K)*ones(4, 1)), 2, 2);
end
end
