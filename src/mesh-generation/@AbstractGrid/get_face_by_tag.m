function indices = get_face_by_tag(obj, marker_name)
indices = find(obj.faces(3, :) == obj.boundary_tags(marker_name));
end
