function a = double(A, options)
%DOUBLE  Assemble the blocks of an IndexedMatrix into one paged array.
%
%   a = double(A, dims = [])
%
%   Rows are offset by the test id and columns by the trial id of each
%   block; blocks with singleton trailing dimensions are broadcast to the
%   largest. dims prescribes the block size of every id; by default it is
%   read from the blocks.
    arguments
        A IndexedMatrix
        options.dims = []
    end
    id_list = cell2mat(A.id);
    if isempty(id_list)
        a = [];
        return
    end

    if isequal(id_list, [0, 0])
        a = double(A.Arrays{:});
        return;
    elseif ismember(-1, id_list)
        id_list = max(id_list + 1, 1);
    end

    A_size = cell2mat(cellfun(@(x) size(x, 1:7), A.Arrays, "UniformOutput", 0));
    if ~isempty(options.dims)
        [row_dims, col_dims] = deal(options.dims(:));
    else
        [row_dims, col_dims] = deal(zeros(max(id_list(:)),1));
        row_dims(id_list(:, 2)) = A_size(:, 1);
        col_dims(id_list(:, 1)) = A_size(:, 2);
    end
    a = zeros([sum(row_dims) sum(col_dims) max(A_size(:, 3:7), [], 1)]);
    
    row_cum_dims = [0; cumsum(row_dims)];
    col_cum_dims = [0; cumsum(col_dims)];
    max_dims = max(A_size(:, 3:7), [], 1);
    for k = 1:length(A.Arrays)
        i = id_list(k, 1);
        i_ = col_cum_dims(i) + (1:col_dims(i));
        j = id_list(k, 2);
        j_ = row_cum_dims(j) + (1:row_dims(j));
        blk = double(A.Arrays{k});
        blk_dims = size(blk, 3:7);
        if ~isequal(blk_dims, max_dims)
            % broadcast blocks with singleton trailing (element/face) dims
            blk = repmat(blk, [1 1 max_dims ./ max(blk_dims, 1)]);
        end
        a(j_, i_, :, :, :, :, :, :) = blk;
    end

end