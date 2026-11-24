function [matching_arrays, non_matching_arrays] = match(A, B)
%MATCH  Pair the blocks of A and B that have the same id.
%
%   [matching_arrays, non_matching_arrays] = match(A, B)
%
%   matching_arrays is a cell of {A block, B block} pairs;
%   non_matching_arrays lists the blocks of A, then of B, that have no
%   partner of the same id in the other matrix.
    arguments
        A IndexedMatrix
        B IndexedMatrix
    end
    a_id = cell2mat(A.id);
    b_id = cell2mat(B.id);
    A_arrays = A.Arrays;
    B_arrays = B.Arrays;

    [is_matching, matching_id] = ismember(a_id, b_id, "rows");

    A_matching_id = find(is_matching);
    B_matching_id = matching_id(matching_id > 0);
    B_non_matching_id = setdiff(1:length(B.Arrays), B_matching_id);

    A_non_matching = A_arrays(~is_matching);
    B_non_matching = B_arrays(B_non_matching_id);
    non_matching_arrays = [A_non_matching; B_non_matching];
    
    gather_match = @(i) {A_arrays{A_matching_id(i)}, B_arrays{B_matching_id(i)}};
    matching_arrays = arrayfun(gather_match, 1:length(A_matching_id), 'UniformOutput', false);
end
