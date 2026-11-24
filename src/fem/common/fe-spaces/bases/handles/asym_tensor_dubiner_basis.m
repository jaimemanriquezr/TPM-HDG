function P = asym_tensor_dubiner_basis(k, r, s)
%ASYM_TENSOR_DUBINER_BASIS  Dubiner basis times the skew matrix [0 1; -1 0].
%
%   P = asym_tensor_dubiner_basis(k, r, s) is tensor_basis with the single
%   linear map [0 1; -1 0].
P = tensor_basis(k, r, s, ... 
                 "PolynomialBasis", @dubiner_basis, ...
                 "LinearMaps", {[0, 1;-1, 0]});
end