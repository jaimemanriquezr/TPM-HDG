function [s, dim, order, id] = read_signature(obj)
%READ_SIGNATURE  Parse the signature string of an FE expression.
%
%   [s, dim, order, id] = read_signature(obj)
%
%   The signature is split at ":". s is the leading name part; dim and
%   order are the trailing numeric fields (-1 when absent); id is the
%   number in the last name field (-1 when it has none).
    signature_values = split(obj.signature, ":");
    [dim, order] = deal(-1);

    if length(signature_values) > 2
        s = signature_values(1:end-2);
        dim = double(signature_values(end-1));
        order = double(signature_values(end));
    elseif length(signature_values) > 1
        s = signature_values(1:end-1);
        dim = double(signature_values(end));
    else
        s = signature_values;
    end
    id = double(extract(s(end), digitsPattern));
    if isempty(id)
        id = -1;
    end
end