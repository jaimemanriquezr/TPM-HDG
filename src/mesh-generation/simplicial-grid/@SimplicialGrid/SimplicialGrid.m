classdef SimplicialGrid < AbstractGrid
    %SimplicialGrid  Conforming triangular mesh of a 2-D domain.
    %
    %   Stores the coordinates of the vertices, the elements (vertex
    %   triples), the faces, the boundaries (face triples of each element)
    %   and the incidence table (element, side, face) of each face, the
    %   orientation and permutation of each element's faces with respect
    %   to the reference element, the normals, the element areas (volume),
    %   the face lengths (area) and the HDG stabilization tau
    %   (add_stabilization).
    %
    %   Boundary faces carry tags in the hidden boundary_tags dictionary
    %   ("Dirichlet", "Neumann", "Interface", "Inflow", ...); dirfaces,
    %   neufaces and intfaces are dependent views of it. Tag faces with
    %   apply_tag, get_face_by_tag or tag_by_curves. Carve a background
    %   grid against an interface with split (or immerse_mesh) and extract
    %   element subsets with submesh or mesh(element_list).
    %
    %   TPM transfer-path data live in the induced partitions of
    %   AbstractCurve.get_partitions, not in the mesh.
    %
    %See also: StructuredSimplicialGrid, split, submesh, apply_tag, fitted_split
    properties (Constant)
        geometric_dimension = 2
        ElementType = "Simplex";
    end

    properties
        coordinates     (2,:) double {mustBeNumeric}    % physical coord.
        elements        (3,:) {mustBeNumeric}           % global (vertex, element)
        faces           (3,:) {mustBeNumeric}           % global (vertex, face)
        boundaries      (3,:) {mustBeNumeric}           % global (face, element)
        incidence       (2, 2, :) {mustBeNumeric}       % global (element, side, face)
        orientation     (3,:) {mustBeNumeric}           % integer -1, 1
        permutation     (3,:) {mustBeNumeric}           % integer 1, 2
        normals         (6,:) double {mustBeNumeric}
        volume          (1,:) double {mustBeNumeric}
        area            (1,:) double {mustBeNumeric}
        tau             (3,:) double {mustBeNumeric}
    end

    methods
        function obj = SimplicialGrid(input_vertices, input_connectivity, options)
            arguments
                input_vertices (2,:) double = []
                input_connectivity (:,:) double = []

                options.filename string = string.empty;
                options.mesh_generator {mustBeKnownMesher} = "mesh2d"; % Triangle, [mesh2d], file, array
                options.element_diameter double = .25;
                options.refine_boundaries = true;

                %%% only for Triangle
                options.triangle_switches string = "Q";
                % q: quality mesh, Q: quiet, Y: do not subdivide boundary edges

                %%% only for mesh2d
                % Spatially varying size target hfun(PP), PP N-by-2, as
                % accepted by refine2; empty = the scalar element_diameter
                options.size_function function_handle = function_handle.empty;
                options.mesh2d_rho2 double = 1.025;
                options.mesh2d_siz1 double = 1.333;
                options.mesh2d_siz2 double = 1.300;
                options.mesh2d_kind string = "delfront"; % [delfront], delaunay
                options.mesh2d_disp = Inf; % Inf = quiet, +10 = standard
            end
            refine_boundaries = options.refine_boundaries;

            if isempty(input_vertices)
                if ~isempty(options.filename)
                    file = fopen(options.filename + ".poly");
                    num = fscanf(file, "%i", 4);
                    v_array = fscanf(file, "%f", [3, num(1)]);
                    v_array(1, :) = [];

                    num = fscanf(file, "%i", 2);
                    c_array = fscanf(file, "%f", [4, num(1)]);
                    c_array(1, :) = [];
                    fclose(file);

                    vertices = v_array;
                    connectivity = c_array;
                else
                    return
                end
            else
                v_1toN = 1:size(input_vertices, 2);
                default_markers = (v_1toN).^0;
                default_numbering = [v_1toN; circshift(v_1toN, -1)];

                vertices = input_vertices;
                connectivity = [default_numbering; default_markers];
            end

            switch size(input_connectivity, 1)
                case 4
                    refine_boundaries = false;
                    h_list = input_connectivity(4, :);
                    h_list(h_list == 0) = options.element_diameter;
                    [vertices, connectivity] = refine_polygon(vertices, ...
                        input_connectivity(1:2, :), input_connectivity(3, :), h_list);
                case 3
                    connectivity = input_connectivity;
                case 2
                    connectivity(1:2, :) = input_connectivity;
                case 1
                    connectivity(3, :) = input_connectivity;
            end

            switch options.mesh_generator
                case "Triangle"
                    if isempty(options.filename)
                        filename = "mesh";
                    else
                        filename = options.filename;
                    end
                    write_polyfile(filename, vertices, connectivity)

                    % the Triangle executable: $TRIANGLE_BIN if set, else
                    % "triangle" on the system PATH
                    tricmd = string(getenv("TRIANGLE_BIN"));
                    if tricmd == ""
                        tricmd = "triangle";
                    end
                    tricmd = tricmd + " -pqe" + options.triangle_switches;
                    if ~refine_boundaries
                        tricmd = tricmd + "Y";
                    end

                    a = get_polygon_area(vertices) * options.element_diameter^2 / 2;
                    [status, message] = system(sprintf(tricmd + "a%.16f %s.poly", a, filename));
                    if status ~= 0
                        error("SimplicialGrid:triangle", ...
                            "Triangle failed (exit status %d; set TRIANGLE_BIN to the " + ...
                            "executable if it is not on the PATH):\n%s", status, message);
                    end

                    %%% reading poly file %%%
                    node = read_poly_block(filename, 4, "node");
                    obj.coordinates = node(1:2, :);
                    obj.elements = read_poly_block(filename, 3, "ele");
                    obj.faces = read_poly_block(filename, 2, "edge");
                case "mesh2d"
                    if refine_boundaries
                        mesh_2d_ref1 = "refine";
                    else
                        mesh_2d_ref1 = "preserve";
                    end
                    if isempty(options.size_function)
                        hfun = options.element_diameter;
                    else
                        hfun = options.size_function;
                    end
                    [node, edge, elem] = refine2(...
                        vertices.', ...
                        connectivity(1:2,:).', ...
                        {}, hfun, ...
                        rho2 = options.mesh2d_rho2, ref1 = mesh_2d_ref1, ...
                        siz1 = options.mesh2d_siz1, siz2 = options.mesh2d_siz2, ...
                        disp = options.mesh2d_disp);
                    edge = tricon2(elem, edge);
                    edge(:, 3:4) = [];

                    boundary_edge = edge(edge(:,3) ~= 0, 1:2);
                    boundary_graph = graph(boundary_edge(:,1), boundary_edge(:,2));
                    for i = 1:size(connectivity, 2)
                        source = connectivity(1, i);
                        target = connectivity(2, i);
                        path = shortestpath(boundary_graph, source, target);
                        subpaths = [path(1:end-1); path(2:end)].';
                        [~, idx1] = ismember(subpaths, edge(:,1:2), "rows");
                        [~, idx2] = ismember(subpaths, fliplr(edge(:,1:2)), "rows");
                        idx = idx1 + idx2;
                        edge(idx, 3) = connectivity(3, i);
                    end
                    obj.coordinates = node.';
                    obj.elements = elem.';
                    obj.faces = edge.';
                case "file"
                    %%% reading poly file %%%
                    filename = options.filename;
                    node = read_poly_block(filename, 4, "node");
                    obj.coordinates = node(1:2, :);
                    obj.elements = read_poly_block(filename, 3, "ele");
                    obj.faces = read_poly_block(filename, 2, "edge");
                case "array"
                    obj.coordinates = input_vertices;
                    obj.elements = input_connectivity;
                    obj.faces = SimplicialGrid.get_faces(input_connectivity);
            end

            [obj.boundaries, ...
                obj.orientation, ...
                obj.incidence, ...
                obj.permutation] = SimplicialGrid.get_boundaries(obj.faces, obj.elements);

            x_K = reshape(obj.coordinates(1, obj.elements), 3, []);
            x_Kp1 = circshift(x_K, -1);
            y_K = reshape(obj.coordinates(2, obj.elements), 3, []);
            y_Kp1 = circshift(y_K, -1);
            obj.volume = .5 * sum(x_K.*y_Kp1 - x_Kp1.*y_K);

            dx_e = diff(reshape(obj.coordinates(1, obj.faces(1:2, :)), 2, []));
            dy_e = diff(reshape(obj.coordinates(2, obj.faces(1:2, :)), 2, []));
            obj.area = sqrt(dx_e.^2 + dy_e.^2);

            n_by_face = [dy_e; -dx_e];
            obj.normals = reshape(n_by_face(:, obj.boundaries), 2*3, []);

            obj.tau = (obj.boundaries).^0;
        end
    end

    methods (Static)
        function out = empty()
            out = SimplicialGrid();
        end

        function faces = get_faces(elements)
            arguments
                elements (3, :)
            end
            edges = permute(cat(3, elements, circshift(elements, 1)), [3 1 2]);
            sorted_edges = reshape(sort(edges), 2, []);

            [unique_edges, ~, edge_indices] = unique(sorted_edges.', 'rows');

            is_interior_edge = (accumarray(edge_indices, 1) == 2);
            unique_edges = unique_edges.';
            bdy_markers = zeros(1, size(unique_edges, 2));

            bdy_markers(~is_interior_edge) = 1;

            faces = [unique_edges; bdy_markers];
        end

        function [boundaries, orientation, incidence, permutation] = get_boundaries(faces, elements)
            [orientation, permutation] = deal(ones(size(elements)));
            orientation = [1; -1; 1] .* orientation;

            % parametrization: e1 = x2->x3, e2 = x1->x3, e3 = x1->x2
            vertex_param = [2 1 1; 3 3 2];
            edge_parametrized = reshape(elements(vertex_param, :), 2, []).';
            edge_faces = faces(1:2, :).';

            % encoding edge into scalars (orientation independent) (x1, x2) -> key(x1,x2)
            n_v  = max([edge_faces(:); elements(:)], [], 'all');
            key_parametrized = min(edge_parametrized, [], 2) + n_v * max(edge_parametrized, [], 2);
            key_faces = min(edge_faces, [], 2) + n_v * max(edge_faces, [], 2);
            [~, location_faces] = ismember(key_parametrized, key_faces);
            assert(all(location_faces > 0), "element edge missing from faces")

            % check reordered edge with original unencoded (orientation dependent)
            edge_is_correct = all(edge_parametrized == edge_faces(location_faces, :), 2);
            orientation(~edge_is_correct) = -1 * orientation(~edge_is_correct);
            permutation(~edge_is_correct) = 2;

            boundaries = reshape(location_faces, 3, []);

            incidence = SimplicialGrid.get_incidence(boundaries);
        end
    end
end

function write_polyfile(name, coords, conn)
arguments
    name    (1,1)  string
    coords  (2,:)  double
    conn    (3,:)  double
end
mesh_dir = fileparts(name);
if ~isfolder(mesh_dir) && (mesh_dir ~= "")
    mkdir(mesh_dir);
end

sc = size(coords,2);
se = size(conn,2);
fID = fopen(strcat(name,".poly"), "w");
fprintf(fID, "%4.1i", [sc 2 0 0]);
fprintf(fID, "\n");
fprintf(fID, "%4.1i %18.16f %18.16f\n", [1:sc; coords]);
fprintf(fID, "\n");
fprintf(fID, "%4.1i", [se 1]);
fprintf(fID, "\n");
fprintf(fID, "%4.1i %4.1i %4.1i %4.1i\n", [1:se; conn]);
fprintf(fID, "\n");
fprintf(fID, "%4.1i",0);
fprintf(fID, "\n");
fprintf(fID, "%4.1i",[1 0 0]);
fclose(fID);
end

function array = read_poly_block(name, len, ext)
file  = fopen(name + ".1." + ext);
num   = fscanf(file, "%i", len);
array = fscanf(file, "%f", [4, num(1)]);
array(1, :) = [];
fclose(file);
end

function area = get_polygon_area(vertices)
x = vertices(1, :);
x_plus_1 = circshift(x, -1);
y = vertices(2, :);
y_plus_1 = circshift(y, -1);
area = reshape(.5 * sum((y + y_plus_1) .* (x - x_plus_1)), [], 1);
end

function mustBeKnownMesher(value)
if ~ismember(value, ["Triangle", "mesh2d", "file", "array"])
    eidType = "mustBeKnownMesher:notKnownMesher";
    msgType = "Name of mesher must be 'Triangle', 'mesh2d', 'file' or 'array'.";
    error(eidType,msgType)
end
end
