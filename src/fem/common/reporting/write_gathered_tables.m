function write_gathered_tables(tables_dir, latex_dir, stem, run_header, ...
    error_by_degree, res_by_degree, nudev_by_degree, mesh_list, degree_list, var_names)
% Gathered tables: every degree completed so far, in one file.
fID = fopen(tables_dir + stem + ".txt", 'w');
fprintf(fID, "%s\n", run_header);
for jj = 1:numel(degree_list)
    fprintf(fID, "\nk = %i\n", degree_list(jj));
    if isempty(nudev_by_degree), nu_page = []; else, nu_page = nudev_by_degree(:, jj); end
    print_error(error_by_degree(:, :, jj), mesh_list, ...
        residue = res_by_degree(:, jj), ...
        normal_deviation = nu_page, names = var_names, fID = fID);
end
fclose(fID);

fID = fopen(latex_dir + stem + ".tex", 'w');
fprintf(fID, "%% %s\n", run_header);
fprintf(fID, "%% k = %s\n", join(string(degree_list), ", "));
print_error_latex(error_by_degree, mesh_list, ...
    names = var_names, degrees = degree_list, ...
    normal_deviation = nudev_by_degree, label_style = "error", fID = fID);
fclose(fID);

% companion table that also prints the normal-deviation angle
if ~isempty(nudev_by_degree)
    fID = fopen(latex_dir + stem + "_with_angles.tex", 'w');
    fprintf(fID, "%% %s\n", run_header);
    fprintf(fID, "%% k = %s\n", join(string(degree_list), ", "));
    print_error_latex(error_by_degree, mesh_list, ...
        names = var_names, degrees = degree_list, ...
        normal_deviation = nudev_by_degree, with_angles = true, ...
        label_style = "error", fID = fID);
    fclose(fID);
end

fprintf("gathered tables: %s{.txt,.tex} (k = %s)\n", ...
    stem, join(string(degree_list), ","));
end
