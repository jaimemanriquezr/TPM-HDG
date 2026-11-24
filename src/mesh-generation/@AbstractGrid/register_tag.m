function [obj, tag_number] = register_tag(obj, tag)
% Adds `tag` to boundary_tags (if new) without tagging any face
arguments
    obj AbstractGrid
    tag string
end
tags_dict = obj.boundary_tags;
if isKey(tags_dict, tag)
    tag_number = tags_dict(tag);
else
    % smallest positive integer not yet used as a tag number
    n_tags = length(values(tags_dict));
    tag_number = min(setdiff(1:n_tags + 1, values(tags_dict)));
    obj.boundary_tags(tag) = tag_number;
end
end
