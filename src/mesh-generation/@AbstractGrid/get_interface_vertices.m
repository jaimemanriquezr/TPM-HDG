function v = get_interface_vertices(obj)
v = obj.coordinates(1:2, obj.faces(1:2, obj.intfaces));
end
