function D = diag(A)
%DIAG  Block-diagonal matrix whose i-th diagonal block is A(:,:,:,i) of double(A).
    A = double(A);
    sz = size(A);
    l = length(sz);
    Z = zeros(size(A, 1:l-1));

    n = sz(end);
    D = cellfun(@(c)Z, cell(n), "UniformOutput", false);
    for i = 1:n
        D{i, i} = A(:, :, :, i);
    end
    D = cell2mat(D);
end