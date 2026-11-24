function P = vector_dubiner_basis(k, r, s)
%VECTOR_DUBINER_BASIS  Dubiner basis times the unit vectors [1; 0] and [0; 1].
%
%   P = vector_dubiner_basis(k, r, s) is tensor_basis with these two
%   linear maps.
P = tensor_basis(k, r, s, ... 
                 "PolynomialBasis", @dubiner_basis, ...
                 "LinearMaps", {[1; 0],[0; 1]});
end