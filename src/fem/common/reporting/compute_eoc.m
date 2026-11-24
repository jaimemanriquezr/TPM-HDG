function [eoc, h_array] = compute_eoc(error_array, mesh_array)
%COMPUTE_EOC  Experimental orders of convergence of an error sequence.
%
%   [eoc, h] = compute_eoc(error_array, mesh_array)
%
%   error_array  n_levels-by-n_vars, one row per refinement level and one
%                column per quantity; consecutive rows are compared.
%   mesh_array   either a double array of mesh sizes h, one per level, or
%                a cell array of SimplicialGrid with one row per level; in
%                the latter case the size of a mesh is max(M.area) (the
%                longest face) and h is the maximum of that over each row.
%
%   eoc  (n_levels-1)-by-n_vars, log(e_{i+1}/e_i) / log(h_{i+1}/h_i);
%        0-by-n_vars when there are fewer than two levels.
%   h_array  the mesh sizes used, one per level.
%
%See also: print_error
    arguments
        error_array double
        mesh_array 
    end
    if isa(mesh_array, "double")
        h_array = mesh_array;
    else
        h_array = cellfun(@(M) max(M.area(:)), mesh_array, "UniformOutput", 0);
        h_array = cell2mat(h_array);
        h_array = max(h_array, [], 2);
    end
    % A single level has no order to observe; the general expression
    % below does not degrade gracefully at one row (empty ./ singleton
    % broadcast error)
    if size(error_array, 1) < 2
        eoc = zeros(0, size(error_array, 2));
        return
    end
    e_log = log(error_array(2:end, :) ./ error_array(1:end-1, :));
    h_log = reshape(log(h_array(2:end) ./ h_array(1:end-1)), [], 1);   % one per row pair
    eoc = e_log ./ h_log;
end
