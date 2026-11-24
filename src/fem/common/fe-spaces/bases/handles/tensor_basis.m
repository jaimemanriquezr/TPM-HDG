function P = tensor_basis(k, r, s, options)
%TENSOR_BASIS  A scalar basis times a list of constant matrices.
%
%   P = tensor_basis(k, r, s, PolynomialBasis = @dubiner_basis, ...
%                    LinearMaps = {eye(2)})
%
%   Each function of the scalar basis is multiplied by every matrix of
%   LinearMaps; the blocks are concatenated map by map along dimension 2.
%   A PolynomialBasis taking more than two arguments is called as
%   (k, r, s), otherwise as (k, r).
    arguments
        k double
        r double
        s double = 0;
        options.PolynomialBasis = @dubiner_basis;
        options.LinearMaps cell = {eye(2)};
    end
    if nargin(options.PolynomialBasis) > 2
        if nargin < 3
            r = r(1, :);
            s = r(2, :);
        end
        p = options.PolynomialBasis(k, r, s);
    else
        p = options.PolynomialBasis(k, r);
    end
    pp = cellfun(@(A) p.*reshape(A, [1 1 size(A)]), ...
                options.LinearMaps, "UniformOutput", false);
    P = cat(2, pp{:});
end