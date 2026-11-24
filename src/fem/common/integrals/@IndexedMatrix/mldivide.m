function x = mldivide(A, b)
%MLDIVIDE  Page-wise solve A\b.
%
%   When b is an IndexedMatrix whose blocks all carry an id of -1 (global
%   right-hand-side blocks), they are first summed into one right-hand
%   side by row offset; otherwise the solve is double(A)\double(b).
    arguments
        A IndexedMatrix
        b
    end
    if isa(b, 'IndexedMatrix') && all(cellfun(@(x) any(x==-1), b.id))
        sz = cellfun(@(x) size(x, 2), A.Arrays);
        A_id = cellfun(@(x) x(1), A.id);
        b_id = cellfun(@(x) x(2), b.id);
        sz_list(A_id) = sz;
        sz_list = [0 sz_list];
        
        double_A = double(A);

        b0 = 0*double_A(:, 1, :, :, :, :);
        for i = 1:length(b.Arrays)
            id = b_id(i);
            double_b = double(b.Arrays{i});

            b0_section = b0(sz_list(id) + (1:sz_list(id+1)), :, :, :, :, :);
            b0(sz_list(id) + (1:sz_list(id+1)), :, :, :, :, :) = b0_section + double_b;
        end
        x = pagemldivide(double_A, b0);
    else
        x = pagemldivide(double(A), double(b));
    end
end
