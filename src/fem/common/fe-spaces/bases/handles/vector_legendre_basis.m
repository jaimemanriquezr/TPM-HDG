function L = vector_legendre_basis(k, s);
%VECTOR_LEGENDRE_BASIS  Legendre basis times the unit vectors [1; 0] and [0; 1].
%
%   L = vector_legendre_basis(k, s) is tensor_basis with these two linear
%   maps.
    L = tensor_basis(k, s, ...
                    "PolynomialBasis", @legendre_basis, ... 
                    "LinearMaps", {[1; 0],[0; 1]});
end