classdef WeakEssentialBoundaryConditions < BoundaryConditions
    %WEAKESSENTIALBOUNDARYCONDITIONS  Dirichlet data imposed weakly on the trace.
    %
    %   bc = WeakEssentialBoundaryConditions(values, tags)
    %
    %   values(i) (Datum) is imposed on the boundary faces tagged tags(i)
    %   through its L2-projection equations onto the face-supported
    %   (trace) global space, assembled as a matrix and a load rather than
    %   by fixing dofs. Vector traces are handled component by component.
    %
    %   [dofs, H, g] = evaluate(bc, problem) returns the global dof indices
    %   of the constrained faces, the sparse matrix H and the load vector g
    %   of the projection equations.
    %
    %See also: BoundaryConditions, StrongEssentialBoundaryConditions
    methods
        function obj = WeakEssentialBoundaryConditions(values, tags)
            arguments
                values Datum = Datum.empty;
                tags string = string.empty;
            end
            obj = obj@BoundaryConditions(values, tags);
        end

        function [dofs, values_H, values_g] = evaluate(obj, problem)
            % essential data attaches to the skeleton-supported (trace)
            % prefix of global_spaces; the system is sized over all supports
            L = global_dof_layout(problem.GlobalSpaces, problem.Mesh);
            n_dof_total = L.n_dof;
            trace_spaces = problem.GlobalSpaces(L.support == "face");
            trace_el = [trace_spaces.Basis];
            n_bc = length(obj.values);
            [weak_bc_dofs, bc_H_values, bc_g_values] = deal(cell(n_bc, 1));
            for ix = 1:n_bc
                bc_faces = get_face_by_tag(problem.Mesh, obj.tags(ix));
                weak_bc_dofs{ix} = get_global_dof(trace_spaces, subset=bc_faces);

                if isa(trace_el, "TensorFEBasis")
                    % vector trace space: assemble component-wise with the
                    % scalar machinery; per-face layout [comp_1; comp_2; ...]
                    scalar_basis = trace_el.ScalarBasis;
                    n_comp = numel(trace_el.LinearMaps);
                    comps = components(obj.values(ix), n_comp);

                    [~, ~, H, g] = project(comps{1}, problem.Mesh, scalar_basis, faces=bc_faces);
                    % the 'ds' facet measure returns one copy per side; keep side 1
                    H = double(H(:, :, :, :, :, 1));
                    g_c = cell(n_comp, 1);
                    g_c{1} = double(g(:, :, :, :, :, 1));
                    for c = 2:n_comp
                        [~, ~, ~, g] = project(comps{c}, problem.Mesh, scalar_basis, faces=bc_faces);
                        g_c{c} = double(g(:, :, :, :, :, 1));
                    end
                    bc_H_values{ix} = pagekron(eye(n_comp), reshape(H, size(H, 1), size(H, 2), []));
                    g_c = cellfun(@(gc) reshape(gc, size(gc, 1), 1, []), g_c, 'UniformOutput', false);
                    bc_g_values{ix} = cat(1, g_c{:});
                else
                    [~, ~, H, g] = project(obj.values(ix), problem.Mesh, trace_el, faces=bc_faces);
                    % the 'ds' facet measure returns one copy per face side
                    % (identical on boundary faces); keep side 1
                    bc_H_values{ix} = double(reshape(H(:, :, :, :, :, 1), size(H, 1), size(H, 2), []));
                    bc_g_values{ix} = double(reshape(g(:, :, :, :, :, 1), size(g, 1), 1, []));
                end
            end
            dofs = [weak_bc_dofs{:}];

            weak_rows = repmat(permute(dofs, [1 3 2]), [1 size(dofs, 1)]);
            weak_cols = permute(weak_rows, [2 1 3]);
            weak_H = cat(3, bc_H_values{:});
            values_H = sparse(weak_rows(:), weak_cols(:), weak_H(:), n_dof_total, n_dof_total);
            values_g = cat(3, bc_g_values{:});
            values_g = values_g(:);
        end

    end
end
