function elements = prune_disconnected_elements(bg_grid, elements)
%-------------------------------------------------------------------------
%Description: Connectivity filter for an active mesh carved out of a
%             background grid. Keeps the face-connected component of
%             largest total area and drops every other component: a
%             connected subdomain must yield a connected mesh.
%
%             Disconnected islands arise when cut elements are assigned
%             to a side (see split's cut_elements) across a stretch of
%             the subdomain thinner than the grid: submesh would tag each
%             island's whole perimeter "Interface", whose projections
%             double-cover Sigma and break the induced partition
%             (get_partitions:foldedProjection).
%
%             Dropped islands are given to neither side, so their region
%             joins the gap band bridged by the transfer paths.
%-------------------------------------------------------------------------
arguments (Input)
    bg_grid   SimplicialGrid  % background grid
    elements  double          % candidate element list
end
arguments (Output)
    elements double  % largest connected component, in the original order
end
elements = elements(:);
if numel(elements) < 2
    return
end

% A global face index appearing twice among the subset's element
% boundaries is a face shared by two of its elements: the dual-graph
% edge list, in local (subset) element indices.
B = bg_grid.boundaries(:, elements);
[sorted_faces, order] = sort(B(:));
shared = find(diff(sorted_faces) == 0);
[~, K_first] = ind2sub(size(B), order(shared));
[~, K_second] = ind2sub(size(B), order(shared + 1));

components = conncomp(graph(K_first, K_second, [], numel(elements)));
if max(components) == 1
    return
end

component_area = accumarray(components(:), abs(bg_grid.volume(elements)));
[~, largest] = max(component_area);
keep = components(:) == largest;
warning("SimplicialGrid:prune_disconnected_elements:droppedIslands", ...
    "Active mesh is disconnected: dropped %d element(s) in %d island(s), " + ...
    "keeping the largest component (%d elements).", ...
    nnz(~keep), max(components) - 1, nnz(keep));
elements = elements(keep);
end
