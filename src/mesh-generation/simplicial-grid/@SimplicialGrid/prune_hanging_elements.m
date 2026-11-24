function elements = prune_hanging_elements(bg_grid, elements)
%-------------------------------------------------------------------------
%Description: Attachment filter for an active mesh carved out of a
%             background grid. Drops every element with two or more
%             carve-exposed faces: a "tooth" hanging off the mesh by a
%             single face.
%
%             Teeth arise on crossed cells with Sigma above the centre
%             vertex: the south triangle survives while its siblings are
%             cut away, turning a staircase boundary into a sawtooth.
%             A tooth face deviates ~42.7 degrees from the interface
%             normal (crossed cell, aspect 1), so it passes the 45-degree
%             default of prune_infeasible_elements, yet it has zero
%             interface measure and its trace dofs carry no transmission
%             condition; the coupled solve diverges when teeth are kept.
%
%             Only faces interior to the background grid and unshared
%             within the kept subset count as carve-exposed; faces on the
%             background grid's own boundary do not. Removing a tooth can
%             expose new ones, so this iterates to a fixed point.
%
%See also: prune_infeasible_elements (the normal-deviation filter),
%          prune_disconnected_elements (the connectivity filter)
%-------------------------------------------------------------------------
arguments (Input)
    bg_grid   SimplicialGrid  % background grid
    elements  double          % candidate element list
end
arguments (Output)
    elements double  % subset attached by at least two faces, original order
end
elements = elements(:);
if numel(elements) < 2
    return
end

% A face of the background grid shared by two of its elements is interior;
% one appearing once lies on the domain boundary. Built once: the grid
% does not change as elements are dropped.
face_multiplicity = accumarray(bg_grid.boundaries(:), 1);

% Each sweep removes at least one element or returns, so numel sweeps can
% exhaust the list and one more observes the fixed point.
for sweep = 0:numel(elements)
    B = bg_grid.boundaries(:, elements);

    % Shared within the subset: a global face index appearing twice among
    % the subset's boundaries, exactly as prune_disconnected_elements
    % builds its dual graph.
    [sorted_faces, order] = sort(B(:));
    twice = find(diff(sorted_faces) == 0);
    shared = false(size(B));
    shared(order(twice)) = true;
    shared(order(twice + 1)) = true;

    % Two or more carve-opened faces. A tooth whose third face is domain
    % boundary (at an inlet or outlet) is deliberately kept: removing it
    % destabilizes the coupled solve.
    opened_by_carve = (face_multiplicity(B) == 2) & ~shared;
    hanging = sum(opened_by_carve, 1) >= 2;
    if ~any(hanging)
        return
    end
    elements(hanging) = [];
end
error("SimplicialGrid:prune_hanging_elements:didNotConverge", ...
    "Hanging-element pruning did not reach a fixed point.");
end
