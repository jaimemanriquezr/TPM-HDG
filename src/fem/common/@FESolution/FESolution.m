classdef FESolution
    % Solution container shared by both solve paths: Unknowns.field /
    % Unknowns.hybrid hold FEField arrays (declarative assembly fills them
    % directly, the low-level path through assign_unknowns).
    properties
        Mesh        AbstractGrid    = SimplicialGrid.empty
        Unknowns    struct          = struct()
    end

    properties (SetAccess = private)
        Problem             FEProblem
        Order               double
    end

    methods
        function obj = FESolution(problem, mesh, k)
            arguments
                problem     FEProblem
                mesh        AbstractGrid    = SimplicialGrid.empty
                k           double          = problem.Order
            end
            obj.Problem = problem;
            obj.Mesh    = mesh;
            obj.Order   = k;
        end

        function obj = assign_unknowns(obj, local_dofs, hybrid_dofs)
            % Wrap low-level dof containers (matrices or per-component
            % cells) as FEFields on the problem's declared bases. Cells
            % stack via cat(1, ...) into the entity layout.
            spaces = obj.Problem.get_discrete_spaces(obj.Order);
            obj.Unknowns = struct( ...
                'field',  fields_from_dofs(spaces.volume,  local_dofs, ...
                                           obj.Problem.FieldsLocal,  obj.Mesh), ...
                'hybrid', fields_from_dofs(spaces.hybrid, hybrid_dofs, ...
                                           obj.Problem.FieldsHybrid, obj.Mesh));
        end

function [L2_field, L2_hybrid] = compute_l2_error(obj, field_exacts, hybrid_exacts)
            % L2 errors of all unknowns against the given exact solutions
            % (cell arrays aligned with Unknowns.field / Unknowns.hybrid;
            % entries are Datum-convertible). Missing entries yield the
            % error against zero. Global scalar unknowns (a ConstantBasis,
            % e.g. a zero-mean multiplier) have no L2 error and yield NaN.
            arguments
                obj           FESolution
                field_exacts  cell = {}
                hybrid_exacts cell = {}
            end
            L2_field  = errors_of(obj.Unknowns.field,  field_exacts,  obj.Mesh);
            L2_hybrid = errors_of(obj.Unknowns.hybrid, hybrid_exacts, obj.Mesh);
        end

        function tiled_layout = plot(obj, fig_handle, options)
            %PLOT  Plot the field unknowns of the solution in one tiled layout.
            %
            %   tiled_layout = plot(sol, fig_handle, components=..., view=...)
            %
            %   fig_handle  figure to draw in (default: a new figure), or an
            %               existing tiledlayout to extend.
            %
            %   Name-value options (defaults):
            %     components  "all" | "tiled" | "first" ("all"), forwarded to
            %                 FEField.plot: every component of the tensor- and
            %                 vector-valued unknowns, each unknown's components
            %                 grouped in a nested layout, or the single
            %                 (1,1) component.
            %     view        3 | 2 (3): 3 keeps the surface view; 2 looks
            %                 straight down with equal axes.
            %
            %   tiled_layout is the layout drawn in, titled with the problem
            %   label and the HDG polynomial degree.
            %
            %See also: FEField/plot
            arguments
                obj         FESolution
                fig_handle  = figure
                % Forwarded to FEField.plot: "all" draws every component
                % of the tensor- and vector-valued unknowns, "tiled"
                % groups each unknown's components in a nested layout,
                % "first" keeps the single-component view.
                options.components (1,1) string ...
                    {mustBeMember(options.components, ["all", "tiled", "first"])} = "all"
                % 3 keeps the surface view; 2 looks straight down and
                % equalises the axes, so the geometry is not stretched to
                % the tile's shape.
                options.view (1,1) double {mustBeMember(options.view, [2, 3])} = 3
            end
            if isa(fig_handle, 'matlab.graphics.layout.TiledChartLayout')
                tiled_layout = fig_handle;
            else
                tiled_layout = tiledlayout(fig_handle, 'flow');
            end

            plottables = obj.Unknowns.field;
            for i = 1:length(plottables)
                plot(plottables(i), obj.Mesh, parent = tiled_layout, ...
                    components = options.components);
            end
            % After the loop, not inside it: an unknown now owns several
            % tiles, and gca() would leave all but its last one flat.
            for ax = findall(tiled_layout, Type = "axes").'
                view(ax, options.view);
                if options.view == 2
                    % equal aspect AND tight limits; the z extent is the
                    % field's, so it is only left alone in the 3-D view
                    axis(ax, "image");
                end
            end

            pde = obj.Problem;
            if isempty(tiled_layout.Title.String)
                title_string = sprintf("\\textbf{%s}", pde.Label);
            else
                title_string = tiled_layout.Title.String + " / " + sprintf("\\textbf{%s}", pde.Label);
            end
            space_letter = "P";
            if obj.Mesh.ElementType == "Quad"
                space_letter = "Q";
            end
            subtitle_string = sprintf("HdG-$\\mathcal{%s}_{%i}$", space_letter, obj.Order);
            title(tiled_layout, title_string, ...
                'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');
            subtitle(tiled_layout, subtitle_string, ...
                'Interpreter', 'latex', 'FontSize', 18);
        end
    end
end

function fields = fields_from_dofs(bases, dofs, metadata, mesh)
% Wrap low-level dof containers as FEFields on the declared bases.
fields = FEField.empty(0, 1);
for i = 1:numel(dofs)
    dof = dofs{i};
    if iscell(dof)
        dof = cat(1, dof{:});
    end
    f = FEField(bases(i), dof);
    f.mesh = mesh;
    if i <= numel(metadata)
        f.name = metadata(i).name;
    end
    fields = [fields; f]; %#ok<AGROW>
end
end

function L2 = errors_of(fields, exacts, mesh)
L2 = zeros(numel(fields), 1);
for i = 1:numel(fields)
    if fields(i).basis.GeometricInputDimension == 0
        % a global scalar (ConstantBasis, e.g. the zero-mean multiplier)
        L2(i) = NaN;  continue
    end
    if i <= numel(exacts) && ~isempty(exacts{i})
        e = exacts{i};
        if ~isa(e, 'Datum'); e = Datum(e); end
        L2(i) = fields(i).compute_l2_error(mesh, e);
    else
        L2(i) = fields(i).compute_l2_error(mesh);
    end
end
end
