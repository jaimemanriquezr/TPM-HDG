classdef StrongEssentialBoundaryConditions < BoundaryConditions
    %STRONGESSENTIALBOUNDARYCONDITIONS  Dirichlet data imposed strongly on the trace.
    %
    %   bc = StrongEssentialBoundaryConditions(values, tags, projection_order=10)
    %
    %   values(i) (Datum) is L2-projected onto the face-supported (trace)
    %   global space on the boundary faces tagged tags(i), and the
    %   resulting dofs are fixed. Vector traces are projected component by
    %   component. projection_order is the quadrature order of the
    %   projections.
    %
    %   [dofs, values] = evaluate(bc, problem) returns, per entry, the
    %   global dof indices and the projected values.
    %
    %See also: BoundaryConditions, WeakEssentialBoundaryConditions
    properties
        projection_order (1,1) double = 10   % quadrature order of the data L2-projections
    end
    methods
        function obj = StrongEssentialBoundaryConditions(values, tags, options)
            arguments
                values Datum = Datum.empty;
                tags string = string.empty;
                % the low-level path projects boundary data with order k
                % (k+1 Gauss points); pass projection_order = k to match it
                % dof-for-dof
                options.projection_order (1,1) double = 10
            end
            obj = obj@BoundaryConditions(values, tags);
            obj.projection_order = options.projection_order;
        end
        
        function [dofs, values] = evaluate(obj, problem)
            % essential data attaches to the skeleton-supported (trace)
            % prefix of global_spaces
            trace_spaces = problem.GlobalSpaces( ...
                [problem.GlobalSpaces.Support] == "face");
            trace_el = [trace_spaces.Basis];
            n_bc = length(obj.values);
            [dofs, values] = deal(cell(n_bc, 1));
            for ix = 1:n_bc
                bc_faces = get_face_by_tag(problem.Mesh, obj.tags(ix));
                if isempty(bc_faces)
                    % no face carries the tag (e.g. a mesh bounded by an
                    % interface only); get_global_dof reads an empty
                    % subset as all faces
                    dofs{ix} = zeros(0, 1);  values{ix} = zeros(0, 1);
                    continue
                end
                dofs{ix} = get_global_dof(trace_spaces, subset=bc_faces);

                if isa(trace_el, "TensorFEBasis")
                    % vector trace: project component-wise with the scalar
                    % machinery (Datum.project's vector path is broken);
                    % per-face dof layout [comp_1; comp_2; ...]
                    scalar_basis = trace_el.ScalarBasis;
                    n_comp = numel(trace_el.LinearMaps);
                    comps = components(obj.values(ix), n_comp);
                    d1 = scalar_basis.Dimension;
                    vals = zeros(n_comp*d1, numel(bc_faces));
                    for c = 1:n_comp
                        proj = project(comps{c}, problem.Mesh, scalar_basis, ...
                            order = obj.projection_order);
                        vals((c-1)*d1 + (1:d1), :) = proj.dof(:, bc_faces);
                    end
                    values{ix} = vals(:);
                else
                    bc_proj = project(obj.values(ix), problem.Mesh, trace_el, ...
                        order = obj.projection_order);
                    values{ix} = reshape(bc_proj.dof(dofs{ix}), [], 1);
                end
            end
        end
    end
end
