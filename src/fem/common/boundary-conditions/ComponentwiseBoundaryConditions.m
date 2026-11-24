classdef ComponentwiseBoundaryConditions < BoundaryConditions
%COMPONENTWISEBOUNDARYCONDITIONS  Per-component Dirichlet/Neumann traces.
%   The declarative counterpart of the low-level path's
%   data.boundary_conditions: a list of boundary_condition(tag, x = ...,
%   y = ..., x_data = ..., y_data = ...) specs, each fixing, on the faces
%   carrying `tag`, every trace component either strongly (dirichlet: the
%   L2 skeleton projection of the datum) or naturally (neumann: the datum
%   is the outward-normal traction component (sigma n)_c, tested against
%   the trace basis and added to the load). Components may mix, e.g.
%   boundary_condition("Symmetry", x = "neumann") fixes u_y and leaves the
%   tangential traction free.
%
%   Projections match the low-level path (L2projSkeleton1d and test_faces
%   at the problem's order); dofs are numbered by the face-supported
%   global space. `values`/`tags` of the base class are unused.
%
%Usage:
%   pb.BoundaryConditions = ComponentwiseBoundaryConditions([ ...
%       boundary_condition("Dirichlet", x_data = g1, y_data = g2);
%       boundary_condition("Outflow",   x = "neumann")]);
%
%See also: boundary_condition, componentwise_boundary_conditions,
%          StrongEssentialBoundaryConditions
    properties
        specs (:,1) struct
    end
    methods
        function obj = ComponentwiseBoundaryConditions(specs)
            arguments
                specs (:,1) struct
            end
            obj = obj@BoundaryConditions(Datum.empty, string.empty);
            obj.specs = specs;
        end

        function faces = natural_faces(obj, mesh)
            % every face with at least one neumann component, i.e. the
            % faces the convective do-nothing flux must skip. An
            % unregistered tag errors here, as in evaluate.
            faces = zeros(0, 1);
            for s = obj.specs.'
                if any(s.type == "neumann")
                    faces = [faces; reshape(get_face_by_tag(mesh, s.tag), [], 1)]; %#ok<AGROW>
                end
            end
            faces = unique(faces);
        end

        function [strong_dofs, strong_vals, natural_dofs, natural_vals] = evaluate(obj, problem)
            mesh = problem.Mesh;  k = problem.Order;
            trace_spaces = problem.GlobalSpaces([problem.GlobalSpaces.Support] == "face");
            assert(isscalar(trace_spaces), "ComponentwiseBoundaryConditions:traceSpaces", ...
                "expected one face-supported global space, found %d.", numel(trace_spaces));
            d1 = trace_spaces.Basis.ScalarBasis.Dimension;      % per component, per face
            % the data is projected at the problem's order, so the trace
            % basis must be of that degree
            assert(d1 == k + 1, "ComponentwiseBoundaryConditions:order", ...
                "trace basis has %d dofs per face but the problem's order is %d.", d1, k);
            [strong_dofs, strong_vals, natural_dofs, natural_vals] = deal(zeros(0, 1));
            for s = obj.specs.'
                if ~isKey(mesh.boundary_tags, s.tag)
                    error("ComponentwiseBoundaryConditions:unknownTag", ...
                        "Boundary tag ""%s"" is not registered on the mesh " + ...
                        "(apply_tag it before solving).", s.tag);
                end
                faces = get_face_by_tag(mesh, s.tag);
                if isempty(faces), continue; end
                % [n_comp*d1, n_faces]: component c of a face occupies rows
                % (c-1)*d1 + (1:d1), as StrongEssentialBoundaryConditions
                % lays its values out
                dofs = get_global_dof(trace_spaces, subset = faces);
                % a spec may set fewer components than the trace has (the
                % rest stay natural with zero data), never more
                assert(numel(s.type) <= size(dofs, 1)/d1, ...
                    "ComponentwiseBoundaryConditions:components", ...
                    "spec ""%s"" sets %d components; the trace has %d.", ...
                    s.tag, numel(s.type), size(dofs, 1)/d1);
                for c = 1:numel(s.type)
                    rows = dofs((c-1)*d1 + (1:d1), :);
                    if s.type(c) == "dirichlet"
                        vals = L2projSkeleton1d(mesh, s.data{c}, k, subset = faces);
                        strong_dofs = [strong_dofs; rows(:)];  %#ok<AGROW>
                        strong_vals = [strong_vals; vals(:)];  %#ok<AGROW>
                    else
                        vals = test_faces(mesh, s.data{c}, k, subset = faces);
                        natural_dofs = [natural_dofs; rows(:)];  %#ok<AGROW>
                        natural_vals = [natural_vals; vals(:)];  %#ok<AGROW>
                    end
                end
            end
        end
    end
end
