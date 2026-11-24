function [mesh_list, list_name] = generate_mesh_list(n_list, interface, gap, options)
%-------------------------------------------------------------------------
%Description: Builds a refinement sequence of split mesh pairs and, by
%             default, saves it for the convergence studies to load.
%
%             Split contract: `split_fun` is called as
%                 [mesh_1, mesh_2] = split_fun(interface, [h_1, h_2], gap)
%             i.e. the two mesh sizes travel together as one vector and
%             the gap has its own slot.
%
%             Anything matching that shape works, so an unfitted family is
%             expressed by closing over a background grid:
%                 split_fun = @(I, h, gap) ...
%                     StructuredSimplicialGrid(1/h(1), 1/h(1), "crossed").split(I)
%
%See also: fitted_split, SimplicialGrid.split
%-------------------------------------------------------------------------
    arguments
        n_list double = 5;
        interface AbstractCurve = FlatCurveHorizontal(.5);
        gap function_handle = @(h) 1*h;

        % Named arguments
        options.h0 double = .2;
        options.i_factor double = 1/2;
        options.split_fun function_handle = @fitted_split

        %%% Options
        options.save_matfile = true;
        options.name_matfile = "./meshes/fitted_split";
    end

    mesh_list = cell(n_list, 2);
    gap_str = func2str(gap);
    for i = 1:n_list
        h1 = options.h0 * (options.i_factor)^(i - 1);
        % Dissimilar sizes are the split_fun's business: a family that
        % wants h_2 ~= h_1 closes over the ratio explicitly.
        h2 = h1;
        d = gap(h1);
        fprintf("Mesh sizes: (%6.5g, %6.5g) --- Gap: %s\n", h1, h2, gap_str)
        fprintf("Generating mesh (%i/%i) ...\n", i, n_list);
        [mesh_list{i, :}] = options.split_fun(interface, [h1, h2], d);
        fprintf("Mesh generated.\n")
    end
    if options.save_matfile
        list_name = options.name_matfile + ".mat";
        list_dir = fileparts(list_name);
        if list_dir ~= "" && ~isfolder(list_dir)
            mkdir(list_dir);
        end
        save(list_name, "mesh_list", "-v7.3");
    end
end
