function normals = get_normal_by_element(obj, element, sideNumber)
% outward normal (scaled by the face length) of side
% sideNumber(i) of element(i), one row per element; normals
% is stored per face orientation, (2*n_sides x n_el)
element = element(:);
sideNumber = sideNumber(:) .* ones(size(element));
rows = 2*(sideNumber - 1) + (1:2);
normals = obj.normals(sub2ind(size(obj.normals), rows, [element, element]));
normals = normals .* obj.orientation(sub2ind(size(obj.orientation), sideNumber, element));
end
