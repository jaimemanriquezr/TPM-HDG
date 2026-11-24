function C = mtimes(A, B)
%MTIMES  Product of an IndexedMatrix with a field, a sign or an array.
%
%   An FEField B multiplies each block by its local dofs, giving an
%   IndexedArray of id [-1 test]; B = 1 or -1 scales each block; an empty
%   operand gives 0 times the other; otherwise the product is pagemtimes of the
%   assembled arrays. The operands are swapped if the first is not an
%   IndexedMatrix.
    arguments
        A
        B
    end
    if ~isa(A, "IndexedMatrix")
        [A, B] = deal(B, A);
    end
    if isa(B, 'FEField')
        v = B.local_dof;
        C = mat_vec_prod(A.Arrays{1}, v);
        for i = 2:length(A.Arrays)
            C = C + mat_vec_prod(A.Arrays{i}, v);
        end
    elseif isequal(B, 1) || isequal(B, -1)
        C = IndexedMatrix(cellfun(@(C) B * C, A.Arrays, "UniformOutput", 0));
    elseif isempty(A)
        C = 0 * B;
    elseif isempty(B)
        C = 0 * A;
    else
        C = pagemtimes(double(A), double(B));
    end
end

function C = mat_vec_prod(A, v)
    % a = N_d2 x N_d1 x 1 x 1 x N_elements x 1 x N_faces_per_element
    % For a SimplicialGrid, N_faces_per_element = 3
    id = A.id;
    a = pagesqueeze(double(A), [2 7]);
    a = reshape(a, size(a, 1), size(a, 2), 1, 1, size(a, 5));
    v = reshape(v, [size(v, 1) 1 1 1 size(v, 3)]);
    C = IndexedArray(pagemtimes(a, v), [-1 id(2)]);
end
