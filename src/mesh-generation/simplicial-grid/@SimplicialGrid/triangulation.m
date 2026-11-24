function T = triangulation(obj)
T = triangulation(obj.elements.', obj.coordinates(1, :).', obj.coordinates(2, :).');
end
