function C = pagesqueeze(A, dims)
%PAGESQUEEZE  Merge two dimensions of an array into one.
%
%   C = pagesqueeze(A, dims)
%
%   dims is a pair of dimensions of A. They are merged into a single
%   dimension placed at the position of the smaller one (the smaller
%   index varies fastest); the other dimensions keep their order. If
%   either entry exceeds ndims(A), A is returned unchanged.
    arguments
        A double
        dims (1, 2) = [1 1];
    end
    if any(dims > length(size(A)))
        C = A;
        return
    end
    dims = sort(dims);
    pre_dims = 1:dims(1)-1;
    post_dims = setdiff(dims(1)+1:length(size(A)), dims(2));
    perm_dims = [pre_dims dims post_dims];

    all_sizes = num2cell(size(A));
    pre_sizes = all_sizes(pre_dims);
    post_sizes = all_sizes(post_dims);
    
    C = reshape(permute(A, perm_dims), pre_sizes{:}, [], post_sizes{:});
end