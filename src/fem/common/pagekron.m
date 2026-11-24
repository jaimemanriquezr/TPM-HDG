function C = pagekron(A, B)
%PAGEKRON  Kronecker product of a matrix with each page of an array.
%
%   C = pagekron(A, B)
%
%   A is m-by-n and B is p-by-q-by-P; C is (m*p)-by-(n*q)-by-P with
%   C(:,:,i) = kron(A, B(:,:,i)). For a 2-D B this is kron(A, B).
    arguments
        A double
        B double
    end

    C = cell2mat(cellfun(@(x) x*B, num2cell(A), 'UniformOutput', false));
end