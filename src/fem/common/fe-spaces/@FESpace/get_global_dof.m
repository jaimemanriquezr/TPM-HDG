function dof_list = get_global_dof(spaces, options)
% Global dof blocks of each space on an entity subset, offset
% by the preceding spaces' dof counts.
arguments
    spaces FESpace
    options.subset double = []
end
n = numel(spaces);
dof_list = cell(n, 1);

num_dof = 0;
for i = 1:n
    assert(~isempty(spaces(i).Mesh), "FESpace:unbound", ...
        "space %d is not mesh-bound: enumeration needs basis x mesh", i);
    i_dof = spaces(i).Basis.get_global_dof(spaces(i).Mesh);
    if ~isempty(options.subset) && ~isa(spaces(i).Basis, "ConstantBasis")
        i_dof = i_dof(:, options.subset);
    end
    dof_list{i} = i_dof + num_dof;
    num_dof = num_dof + length(unique(i_dof));
end

if isscalar(dof_list)
    dof_list = dof_list{:};
else
    sz = cellfun(@(x) size(x, 3:length(size(x))), dof_list, 'UniformOutput', false);
    if isequal(sz{:})
        dof_list = repmat(cell2mat(dof_list), [1 n]);
    end
end
end
