function new_obj = mtimes(obj, A)
arguments
    obj
    A
end
if ~isa(obj, "SimplicialGrid")
    [obj, A] = deal(A, obj);
end
new_obj = obj;
new_obj.coordinates = pagemtimes(A, obj.coordinates);
new_obj.volume = abs(det(A)) * obj.volume;

% the normals carry the face length, so they transform with
% the cofactor matrix det(A) A^{-T}
nu = reshape(obj.normals, 2, 3, []);
new_nu = det(A) * pagemldivide(A, "transpose", nu);

% inefficient loop
new_obj.area = obj.area;
for f = 1:size(obj.faces, 2)
    nu_f = new_nu(:, obj.incidence(1, 2, f), obj.incidence(1, 1, f));
    new_obj.area(f) = sqrt(sum(nu_f.^2));
end
new_obj.normals = reshape(new_nu, 3*2, []);
end
