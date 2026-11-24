function new_obj = submesh(obj, element_list)
arguments
    obj SimplicialGrid
    element_list (1, :) double = 1:size(obj.elements, 2);
end

new_obj = obj; % copy

new_obj.tau = obj.tau(:, element_list);
new_obj.orientation = obj.orientation(:, element_list);
new_obj.permutation = obj.permutation(:, element_list);
new_obj.normals = obj.normals(:, element_list);
new_obj.volume = obj.volume(element_list);

active_vertices = obj.elements(:, element_list);
vertices_list = unique(active_vertices);

active_boundaries = obj.boundaries(:, element_list);
faces_list = unique(active_boundaries);

active_faces = obj.faces(1:2, faces_list);
active_bmarker = obj.faces(3, faces_list);
new_obj.area = obj.area(faces_list);

new_vertices = (1:numel(vertices_list)).';
[~, element_loc] = ismember(active_vertices, vertices_list);
[~, face_loc] = ismember(active_faces, vertices_list);

new_faces = (1:numel(active_faces)).';
[~, boundary_loc] = ismember(active_boundaries, faces_list);

new_obj.coordinates = obj.coordinates(:, vertices_list);
new_obj.elements = new_vertices(element_loc);
new_obj.faces = [new_vertices(face_loc); active_bmarker];
new_obj.boundaries = new_faces(boundary_loc);
% the face incidence table refers to old element/face numbers;
% rebuild it for the submesh
new_obj.incidence = SimplicialGrid.get_incidence(new_obj.boundaries);

% retag the faces the cut exposes: boundary faces of the
% submesh that were interior in the parent
boundary_faces = (accumarray(new_obj.boundaries(:), 1) == 1);
parent_count = accumarray(obj.boundaries(:), 1, [size(obj.faces, 2), 1]);
new_inter_faces = boundary_faces & (parent_count(faces_list) == 2);

new_obj.faces(3, new_inter_faces) = obj.boundary_tags("Interface");
end
