function write_field(fID, name, dims, type, value)
%Description: Writes a named field to an open text file in the txt2mat format.
    fprintf(fID, "%s %i %i %s\n", name, dims(:), type);
    switch type
        case "function_handle"
            fprintf(fID, "%s\n", string(value));
        otherwise
            fprintf(fID, "%g\n", value);
    end
end
