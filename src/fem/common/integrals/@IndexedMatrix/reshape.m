function B = reshape(A, varargin)
%RESHAPE  Reshape every block, keeping its id; the blocks must have equal sizes.
    sz = cellfun(@size, A.Arrays, "UniformOutput", false);
    assert(isscalar(sz) || isequal(sz{:}), "Can only reshape homogeneous IndexedMatrix.");
    
    B = cell(size(A.Arrays));
    for i = 1:length(A.Arrays)
        % reshape of the double subclass drops the id: rebuild the block
        B{i} = IndexedArray(reshape(double(A.Arrays{i}), varargin{:}), A.Arrays{i}.id);
    end
    B = IndexedMatrix(B);
end