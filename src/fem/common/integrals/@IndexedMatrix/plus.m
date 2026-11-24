function A_plus_B = plus(A, B)
%PLUS  Sum of two IndexedMatrix objects.
%
%   Blocks of equal id are added, the others are concatenated. A matrix of
%   id [-1 -1] or an empty one is neutral; if either has id [0 0] the sum
%   is that of the assembled arrays.
    arguments
        A IndexedMatrix
        B IndexedMatrix
    end
    if isequal(A.id, {[-1 -1]})
        A_plus_B = B;
        return
    elseif isequal(B.id, {[-1 -1]})
        A_plus_B = A;
        return
    elseif isempty(A)
        A_plus_B = B;
        return
    elseif isempty(B)
        A_plus_B = A;
        return
    elseif isequal(A.id, {[0 0]}) || isequal(B.id, {[0 0]})
        A_plus_B = double(A) + double(B);
        return
    end

    [matching_arrays, non_matching_arrays] = match(A, B);
    matched_arrays = arrayfun(@(i) plus(matching_arrays{i}{:}), 1:length(matching_arrays), 'UniformOutput', false);
    A_plus_B = IndexedMatrix([matched_arrays(:); non_matching_arrays(:)]);
end
