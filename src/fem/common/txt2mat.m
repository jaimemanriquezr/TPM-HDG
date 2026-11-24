function output = txt2mat(filename)
%TXT2MAT  Convert a text file of named arrays and expressions to a MAT file.
%
%   txt2mat(filename)
%   data = txt2mat(filename)
%
%   The file holds records of a header line "name rows cols type" followed
%   by rows*cols lines. For type "double" each line is a number, filled
%   column-major into a rows-by-cols array. For type "function_handle"
%   each line is an expression in x and y, turned into an elementwise
%   handle @(x, y) (the operators * ^ / become .* .^ ./) and collected in
%   a rows-by-cols cell, or a single handle when rows*cols is 1.
%
%   The variables are saved to filename with a trailing ".txt" replaced
%   by ".mat" (the workspace is saved as is, including filename). With an
%   output, the struct loaded from that file is returned.
op_syms = ["*", "^", "/"];
fID = fopen(filename,'r');
while ~feof(fID)
    fields = split(fgetl(fID));
    [name, sz1, sz2, type] = deal(fields{:});
    sz1 = str2double(sz1); sz2 = str2double(sz2);

    str = "";
    switch type
        case "double"
            for n = 1:(sz1*sz2)
                str = str + fgetl(fID) + ";";
            end
            eval_line = sprintf(name + "=zeros(%i,%i); ", sz1, sz2) + ...
                name + "(:) = [" + str + "]";
        case "function_handle"
            for n = 1:(sz1*sz2)
                str = str + "@(x,y) (x.*y).^0 .*(" + replace(string(fgetl(fID)),op_syms,"." + op_syms) + ");";
            end
            eval_line = sprintf(name + "=cell(%i,%i); ", sz1, sz2) + ...
                name + "(:) = {" + str + "}";

            if sz1*sz2 == 1
                eval_line = eval_line + "; " + name + " = " + name + "{1}";
            end
    end
    eval(eval_line + ";");
end
fclose(fID);
clear ans eval_line fID fields n name str sz1 sz2 type op_syms
filename = regexprep(filename, '\.txt$', '');
save(filename + ".mat");
if nargout > 0
    output = load(filename + ".mat");
end
end
