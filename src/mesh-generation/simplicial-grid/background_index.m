function index = background_index(carved, background)
%-------------------------------------------------------------------------
%Description: Which element of `background` each element of `carved` is.
%
%             `split` and `submesh` select a subset of a background grid's
%             elements and renumber everything, keeping no back-pointer, so
%             the correspondence has to be recovered. It is exact rather
%             than approximate: a carved element is a background element,
%             not merely contained in one, so locating its centroid
%             identifies it uniquely.
%
%             Needed when an unfitted hierarchy is built by re-carving a
%             refined background: the parent maps live on the background
%             grids (where uniform_refine produced them), and this is what
%             carries them over to the carved meshes.
%
%Usage:       idx = background_index(fluid_mesh, bg_grid)
%
%See also: SimplicialGrid.split, SimplicialGrid.uniform_refine
%-------------------------------------------------------------------------
arguments
    carved     SimplicialGrid
    background SimplicialGrid
end

centroids = squeeze(mean(reshape( ...
    carved.coordinates(:, carved.elements), 2, 3, []), 2));
index = pointLocation(triangulation(background), centroids.').';

assert(~any(isnan(index)), "background_index:notFound", ...
    "%d carved centroid(s) lie in no background element; the two meshes " + ...
    "are not related by a selection.", sum(isnan(index)));

% A selection is injective. If it is not, the caller has handed us meshes
% from different grids that happen to overlap.
assert(numel(unique(index)) == numel(index), "background_index:notInjective", ...
    "Two carved elements map to the same background element.");

% The areas must agree element by element, which is the check that the
% carved element really is the background element and not just inside it.
drift = max(abs(abs(carved.volume(:)) - abs(background.volume(index).')));
assert(drift < 1e-10, "background_index:areaMismatch", ...
    "Carved and background element areas differ by %.3e; the carved mesh " + ...
    "is not a selection from this grid.", drift);
end
