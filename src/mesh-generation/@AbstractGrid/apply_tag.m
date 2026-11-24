function obj = apply_tag(obj, tag, tag_indicator)
arguments
    obj AbstractGrid
    tag string
    tag_indicator function_handle % boolean output
end
[obj, tag_number] = register_tag(obj, tag);

x = obj.coordinates;
f_list = obj.faces(1:end-1, :);
is_on_tag = all(reshape(tag_indicator(x(:, f_list)), 2, []));

obj.faces(3, is_on_tag) = tag_number;
end
