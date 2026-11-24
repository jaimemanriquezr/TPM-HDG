function print_error_latex(error_list, mesh_list, options)
% LaTeX tabular of errors and observed orders of convergence.
%
% error_list is n_levels x n_vars for a single polynomial degree. Passing
% `degrees` switches to the gathered layout: error_list is then
% n_levels x n_vars x n_degrees and every degree becomes a \multirow block
% of one table, sharing the header.
%
% normal_deviation, when given, adds || nu_h - nu ||_{L^inf} between the
% discrete and exact interface normals as a column after $h$ with its own
% eoc, as in print_error. n_levels x 1 for a single degree, n_levels x
% n_degrees for the gathered layout.
arguments
    error_list double
    mesh_list
    options.names   = string.empty
    options.fID     = 1
    options.degrees double = []
    options.normal_deviation double = []
    % Append the angle theta = 2*asin(|nu_h - nu|/2), in degrees, beside
    % the deviation pair.
    options.with_angles (1,1) logical = false
    % Column-head style. "plain" writes $u$; "error" writes $e(\bs{u})$
    % (\hat{u} becomes \hat{e}(\bs{u})) and expects \bs from the including
    % document's preamble.
    options.label_style (1,1) string ...
        {mustBeMember(options.label_style, ["plain", "error"])} = "plain"
    % Header and per-block labels of the \multirow groups (default: the
    % polynomial degree). Written verbatim, with their own math delimiters.
    options.group_header (1,1) string = "$k$"
    options.group_labels (1,:) string = string.empty
end
n_vars    = size(error_list, 2);
n_degrees = size(error_list, 3);
gathered  = ~isempty(options.degrees);
if gathered && numel(options.degrees) ~= n_degrees
    error("print_error_latex:degreeMismatch", ...
        "degrees has %d entries but error_list has %d pages.", ...
        numel(options.degrees), n_degrees);
end

nu_dev  = options.normal_deviation;
with_nu = ~isempty(nu_dev);
% One column per degree page: the deviation is sampled with a k-dependent
% quadrature, so each \multirow block has its own.
if with_nu && size(nu_dev, 2) ~= n_degrees
    error("print_error_latex:normalDeviationMismatch", ...
        "normal_deviation has %d column(s) but error_list has %d page(s).", ...
        size(nu_dev, 2), n_degrees);
end

if isempty(options.names)
    names = string(arrayfun(@(j) "x_{" + j + "}", 1:n_vars, UniformOutput=false));
else
    names = options.names;
end

col_spec = "c|" + join(repmat("rr", n_vars, 1), "|");
header   = "$h$";
if with_nu
    % Not "rr" + repmat("r", 1, with_angles): repmat to 0 copies gives an
    % empty string array, which empties the whole spec.
    if options.with_angles, nu_cols = "rrr"; else, nu_cols = "rr"; end
    col_spec = "c|" + nu_cols + "|" + join(repmat("rr", n_vars, 1), "|");
    header   = header + " & " + error_label("\nu", options.label_style) + " & eoc";
    if options.with_angles
        header = header + " & $\angle(\nu_h,\nu)$";
    end
end
if gathered
    col_spec = "c|" + col_spec;
    header   = options.group_header + " & " + header;
end
% One label per block; the default is the degree.
group_labels = options.group_labels;
if gathered && isempty(group_labels)
    group_labels = compose("%d", options.degrees);
end
if gathered && numel(group_labels) ~= n_degrees
    error("print_error_latex:groupLabelMismatch", ...
        "group_labels has %d entries but error_list has %d pages.", ...
        numel(group_labels), n_degrees);
end
for j = 1:n_vars
    header = header + " & " + error_label(names(j), options.label_style) + " & eoc";
end
fprintf(options.fID, "\\begin{tabular}{%s}\n\\hline\n", col_spec);
fprintf(options.fID, "%s \\\\\n\\hline\n", header);

for d = 1:n_degrees
    errors = error_list(:, :, d);
    [eoc, h] = compute_eoc(errors, mesh_list);
    eoc = [nan(1, n_vars); eoc]; %#ok<AGROW>
    if with_nu
        % Against the same h as the errors. NaN (printed as ---) on a
        % fitted family, where the deviation is 0 and the ratio is 0/0.
        nu_col = nu_dev(:, d);
        nu_eoc = [nan; compute_eoc(nu_col, h)];
    end
    n_rows = size(errors, 1);
    for i = 1:n_rows
        if ~gathered
            row = "";
        elseif i == 1
            row = sprintf("\\multirow{%d}{*}{%s} & ", n_rows, group_labels(d));
        else
            row = " & ";
        end
        row = row + sprintf("%.2e", h(i));
        if with_nu
            if isnan(nu_eoc(i))
                row = row + sprintf(" & %.2e & ---", nu_col(i));
            else
                row = row + sprintf(" & %.2e & %.2f", nu_col(i), nu_eoc(i));
            end
            if options.with_angles
                % The exact inverse of the norm: |nu_h - nu| = 2 sin(theta/2).
                row = row + sprintf(" & %.1f", 2 * asind(nu_col(i) / 2));
            end
        end
        for j = 1:n_vars
            if isnan(eoc(i, j))
                row = row + sprintf(" & %.2e & ---", errors(i, j));
            else
                row = row + sprintf(" & %.2e & %.2f", errors(i, j), eoc(i, j));
            end
        end
        fprintf(options.fID, "%s \\\\\n", row);
    end
    if gathered, fprintf(options.fID, "\\hline\n"); end
end

if ~gathered, fprintf(options.fID, "\\hline\n"); end
fprintf(options.fID, "\\end{tabular}\n");
end

function label = error_label(name, style)
% Column head for `name`. "plain" is the variable itself; "error" moves a
% hat to the error: \hat{u} becomes \hat{e}(\bs{u}). \bs must be defined
% by the including document.
if style == "plain"
    label = "$" + name + "$";
elseif startsWith(name, "\hat{") && endsWith(name, "}")
    inner = extractBetween(name, strlength("\hat{") + 1, strlength(name) - 1);
    label = "$\hat{e}(\bs{" + inner + "})$";
else
    label = "$e(\bs{" + name + "})$";
end
end
