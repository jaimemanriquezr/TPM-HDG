function [F, x0] = get_affine_map(obj, K, options)
%GET_AFFINE_MAP  Affine maps from the reference square to grid elements.
%
%   [F, x0] = get_affine_map(mesh, K, inverse = false, output = "array")
%
%   Element K (default all elements) is the image of [0, 1]^2 under
%   x = B*[r; s] + x0, with B = diag(hx, hy) and x0 its south-west
%   corner. F is B (or inv(B) when inverse is true), returned as a
%   2-by-2-by-numel(K) array for output = "array" or as a 2-by-2 cell of
%   column vectors for output = "cell", as in SimplicialGrid; x0 is
%   2-by-numel(K). The empty mesh gives the identity map.
%
%See also: SimplicialGrid/get_affine_map, get_2d_local
arguments
    obj CartesianGrid
    K (:,1) = 1:size(obj.elements, 2)
    options.inverse = false
    options.output = "array"
end

if isempty(obj)
    F = eye(2);
    x0 = [0.; 0.];
    return
end

[i, j] = ind2sub([obj.Nx, obj.Ny], K);
B_11 = obj.hx(i);
B_22 = obj.hy(j);
x0 = [obj.x(i); obj.y(j)];
if options.inverse
    [B_11, B_22] = deal(1 ./ B_11, 1 ./ B_22);
end
z = zeros(size(B_11));
F = [B_11(:), z(:), z(:), B_22(:)];

switch options.output
    case "array"
        F = reshape(F.', 2, 2, []);
    case "cell"
        F = reshape(mat2cell(F(:), length(K)*ones(4, 1)), 2, 2);
end
end
