function print_error(error_list, mesh_list, options)
%PRINT_ERROR  Print an error table with EOC columns.
%
%   print_error(error_list, mesh_list, Name = Value, ...)
%
%   error_list  n_levels-by-n_vars errors, one row per level (as in
%               compute_eoc).
%   mesh_list   mesh sizes or cell array of SimplicialGrid (compute_eoc's
%               mesh_array); the first column of the table is h.
%
%   Options:
%     names            column labels, one per quantity.
%     format_error     printf format of h and of the error columns.
%     format_eoc       printf format of the EOC columns.
%     residue          n_levels-by-1 values printed as an extra column of
%                      their own, left of the error columns, in
%                      format_error (no EOC column).
%     normal_deviation n_levels values of the interface normal deviation
%                      max-norm, printed as a column "nu_h-nu" (with its
%                      EOC if print_eoc) beside the errors.
%     print_eoc        logical; false omits the EOC columns.
%     fID              file identifier the table is written to.
%
%See also: compute_eoc
    arguments
    error_list double
    mesh_list
    options.names = string.empty
    options.format_error = "%9.3e";
    options.format_eoc = "%4.1f";
    
    options.residue = [];
    % || nu_h - nu ||_{L^inf} between each mesh's discrete interface-face
    % normal and the exact interface normal at the patch quadrature
    % points, printed beside the errors as a column of its own.
    options.normal_deviation = [];
    options.print_eoc = true;
    
    options.fID = 1;
    end
    [n_rows, n_vars] = size(error_list, 1:2);
    [eoc, h] = compute_eoc(error_list, mesh_list);

    if isempty(options.names)
        name_list = repmat("", [n_vars, 1]);
    else
        name_list = options.names;
    end
    p_array = error_list;

    format_list = repmat(options.format_error, [1 n_vars]);
    if options.print_eoc
        name_list = [name_list.'; repmat("eoc", 1, length(name_list))];
        format_list = [format_list; repmat(options.format_eoc, [1 n_vars])];
        
        % nan(1, n_vars), not nan*eoc(1, :): with a single row compute_eoc
        % returns 0-by-n_vars and indexing its first row errors.
        p_array = cat(3, error_list, [nan(1, n_vars); eoc]);
        p_array = permute(p_array, [1 3 2]);
        p_array = reshape(p_array, size(p_array, 1), [], 1);
    end

    if ~isempty(options.residue)
        name_list = ["residue"; name_list(:)];
        format_list = string([options.format_error; format_list(:)]);
        p_array = [options.residue, p_array];
    end
    if ~isempty(options.normal_deviation)
        nu_dev     = options.normal_deviation(:);
        nu_names   = "nu_h-nu";
        nu_formats = options.format_error;
        nu_columns = nu_dev;
        if options.print_eoc
            % Against the same h as the errors. On a fitted family
            % nu_h - nu is 0 by construction (0/0), so the column is NaN.
            nu_names   = [nu_names; "eoc"];
            nu_formats = [nu_formats; options.format_eoc];
            nu_columns = [nu_dev, [nan; compute_eoc(nu_dev, h)]];
        end
        name_list = [nu_names(:); name_list(:)];
        format_list = string([nu_formats(:); format_list(:)]);
        p_array = [nu_columns, p_array];
    end
    name_list = ["h"; name_list(:)];
    format_list = string([options.format_error; format_list(:)]);
    p_array = [h, p_array];
    

    f_e = double(extract(extract(format_list, "%" + digitsPattern), digitsPattern));
    print_line_strings = @(list) print_line(list, "%" + f_e + "s", options.fID);
    print_line_double = @(list) print_line(list, format_list, options.fID);
    hrule = @() print_horizontal_rule(f_e, options.fID);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hrule();
    print_line_strings(name_list);
    hrule();
    for i = 1:n_rows
        print_line_double(p_array(i, :));
    end
    hrule();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function print_line(A, s_format, fID)
    arguments
        A 
        s_format
        fID = 1
    end
    fprintf(fID, '|'); 
    fprintf(fID, join(s_format + "|"), A); 
    fprintf(fID, '\n');
end

function print_horizontal_rule(f_e, fID)
    arguments
        f_e
        fID = 1
    end
    for i = 1:length(f_e)
        fprintf(fID, repmat('-', 1, f_e(i) + 1));
        fprintf(fID, '-');
    end
    fprintf(fID, '\n');
end