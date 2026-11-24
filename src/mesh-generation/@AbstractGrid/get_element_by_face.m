function [element, side_number] = get_element_by_face(obj, face)
arguments
    obj AbstractGrid
    face = 1:size(obj.faces, 2);
end
K = permute(obj.incidence(:, 1, face), [1 3 2]);   % 2 x numel(face)
s = permute(obj.incidence(:, 2, face), [1 3 2]);

% the incidence table duplicates boundary faces; return a single
% entry for them
is_single = K(1, :) == K(2, :) & s(1, :) == s(2, :);

if ~any(is_single)                       % all interior: 2 x n arrays
    element     = K;
    side_number = s;
elseif all(is_single)                    % all boundary: 1 x n arrays
    element     = K(1, :);
    side_number = s(1, :);
else                                     % mixed: cells
    [element, side_number] = deal(cell(size(face)));
    for e = 1:numel(face)
        element{e}     = K(1:2 - is_single(e), e);
        side_number{e} = s(1:2 - is_single(e), e);
    end
end
end
