function [rows, cols, rows_rhs] = get_assembly_dof(spaces, mesh)
% Assembly scatter indices over a (possibly heterogeneous)
% space array: per-space blocks from each element's own local
% enumeration, stitched with running dof offsets.
arguments
    spaces FESpace
    mesh AbstractGrid
end
n = numel(spaces);
rows = cell(n, 1);
cols = cell(1, n);

num_dof = 0;
for i = 1:n
    % dof enumeration is contiguous per face for scalar and
    % vector (tensor) trace elements alike: vector dofs are
    % stored per face as [comp_1; comp_2; ...], matching the
    % pagekron block structure of the assembled matrices
    [r, c] = spaces(i).Basis.get_assembly_dof(mesh);
    rows{i} = r + num_dof;
    cols{i} = c + num_dof;
    num_dof = num_dof + length(unique(r));
end

if isscalar(rows)
    rows = rows{:};
    rows_rhs = rows(:, 1, :, :, 1);
else
    sz = cellfun(@(x) size(x, 3:length(size(x))), rows, 'UniformOutput', false);
    if isequal(sz{:})
        rows = repmat(cell2mat(rows), [1 n]);
        rows_rhs = rows(:, 1, :, :, 1);
    end
end

if isscalar(cols)
    cols = cols{:};
else
    sz = cellfun(@(x) size(x, 3:length(size(x))), cols, 'UniformOutput', false);
    if isequal(sz{:})
        cols = repmat(cell2mat(cols), [n 1]);
    end
end
end

