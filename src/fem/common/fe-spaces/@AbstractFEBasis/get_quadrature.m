function [x_local, w] = get_quadrature(bases, dim)
%GET_QUADRATURE  Reference quadrature rule for a set of bases.
%
%   [x_local, w] = get_quadrature(bases)
%   [x_local, w] = get_quadrature(bases, dim)
%
%   Returns the points x_local and weights w of the reference quadrature
%   rule of degree k + 1, where k is the largest Order in bases (the rule
%   comes from get_2d_reference_quadrature / get_1d_reference_quadrature
%   for the ElementType of the element-supported bases, or of all bases if
%   there are none). dim = 2 gives the element rule, dim = 1 the
%   facet rule; it defaults to the smallest nonzero GeometricInputDimension
%   among the bases (2 if all are scalar-supported).
arguments
    bases AbstractFEBasis
    % scalar-supported bases (gid = 0) carry no geometry: they follow the
    % quadrature of the lowest-dimensional geometric basis in the set
    dim double = min_geometric_dimension(bases);
end
k = max([bases.Order]);
el_type = element_type(bases);
switch dim
    case 2
        [x_local, w] = get_2d_reference_quadrature(k + 1, el_type);
    case 1
        [x_local, w] = get_1d_reference_quadrature(k + 1, el_type);
end
end

function el_type = element_type(bases)
% the element type of the element-supported bases; face and scalar bases
% (LegendreBasis, ConstantBasis) fit any element, so they decide only
% when no element-supported basis is present
on_elements = bases([bases.GeometricInputDimension] == 2);
if isempty(on_elements)
    on_elements = bases;
end
el_type = unique([on_elements.ElementType]);
end

function dim = min_geometric_dimension(bases)
dims = [bases.GeometricInputDimension];
dims = dims(dims > 0);
if isempty(dims)
    dim = 2;   % a set of only scalar-supported bases: element quadrature
else
    dim = min(dims);
end
end
