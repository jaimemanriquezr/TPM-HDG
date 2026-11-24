classdef (Abstract) AbstractGrid < matlab.mixin.indexing.RedefinesParen
    %AbstractGrid  Common parent of the 2-D meshes (SimplicialGrid, CartesianGrid).
    %
    %   Subclasses must define the mesh tables as properties, with their
    %   own size validators (which MATLAB does not allow on inherited
    %   abstract properties, so they are not declared here):
    %     coordinates   (2, n_vertices)
    %     elements      (vertex, element)
    %     faces         (vertex 1, vertex 2, tag; face)
    %     boundaries    (face, element): the face of each element side
    %     incidence     (neighbour, [element side], face)
    %     orientation   (side, element), integer -1, 1: the outward normal
    %                   of the side is orientation .* the face normal
    %     permutation   (side, element), integer 1, 2: 2 when the side
    %                   runs against the face parametrization
    %     normals       (2*n_sides, element): the face normal
    %                   [dy; -dx] of the face of each side, v1 -> v2,
    %                   scaled by the face length
    %     volume        element areas;  area   face lengths
    %     tau           HDG stabilization (side, element)
    %   and implement:
    %     empty()                  the empty mesh (from RedefinesParen)
    %     get_faces(elements)      faces (vertex 1, vertex 2, tag) of a
    %                              connectivity, boundary faces tagged 1
    %     get_boundaries(faces, elements)   [boundaries, orientation,
    %                              incidence, permutation], the faces of
    %                              each element side and their orientation
    %     get_affine_map(obj, K, ...)       x = F*[r; s] + x0 per element
    %     get_2d_local(obj, xy, elements)   inverse of the element maps
    %     get_2d_local_boundary(obj, t, elements, sides)   reference
    %                              points of face parameters
    %     submesh(obj, element_list)        used by mesh(element_list)
    %
    %   Inherited: boundary tags (get_face_by_tag, register_tag, apply_tag,
    %   dirfaces, neufaces, intfaces, get_interface_vertices), incidence
    %   queries (get_element_by_face), the geometry built on the affine
    %   maps (get_2d_physical, get_physical_points, get_normal_by_element,
    %   exterior_unit_normals), add_stabilization, plot, vertices, size,
    %   and mesh(element_list) indexing.
    %
    %See also: SimplicialGrid, CartesianGrid
    properties (Constant, Abstract)
        ElementType
        geometric_dimension
    end

    properties (Hidden)
        boundary_tags = dictionary(["Dirichlet", "Neumann", "Interface", "Inflow"], 1:4)
    end

    properties (Hidden, Dependent)
        dirfaces        (1,:) {mustBeNumeric}           % global numbering
        neufaces        (1,:) {mustBeNumeric}           % global numbering
        intfaces        (1,:) {mustBeNumeric}           % global numbering
        vertices        (2, :, :) {mustBeNumeric}       % coords by ele
        exterior_unit_normals (2, 1, :, :) {mustBeNumeric}  % (2, 1, side, element)
    end

    methods (Abstract)
        new_obj = submesh(obj, element_list)
        [F, x0] = get_affine_map(obj, K, options)
        [rs_points, element_list] = get_2d_local(obj, xy_points, element_list)
        rs_points = get_2d_local_boundary(obj, t_points, element_list, side_number)
    end

    methods (Abstract, Static)
        faces = get_faces(elements)
        [boundaries, orientation, incidence] = get_boundaries(faces, elements)
    end

    methods
        indices = get_face_by_tag(obj, marker_name)
        [obj, tag_number] = register_tag(obj, tag)
        obj = apply_tag(obj, tag, tag_indicator)
        v = get_interface_vertices(obj)
        [element, side_number] = get_element_by_face(obj, face)
        obj = add_stabilization(obj, tau)
        p = plot(obj, y, options)
        xy_points = get_2d_physical(obj, rs_points, element_list)
        x_physical = get_physical_points(obj, x_local)
        normals = get_normal_by_element(obj, element, sideNumber)

        function N = get.exterior_unit_normals(obj)
            scaled_orientation = obj.orientation ./ (obj.area(obj.boundaries));
            N = reshape(obj.normals .* kron(scaled_orientation, [1; 1]), ...
                2, 1, size(obj.elements, 1), []);
        end

        function value = get.dirfaces(obj)
            value = get_face_by_tag(obj, "Dirichlet");
        end

        function value = get.neufaces(obj)
            value = get_face_by_tag(obj, "Neumann");
        end

        function value = get.intfaces(obj)
            value = get_face_by_tag(obj, "Interface");
        end

        function value = get.vertices(obj)
            value = reshape(obj.coordinates(:, obj.elements), 2, size(obj.elements, 1), []);
        end
    end

    methods (Static, Access = protected)
        function incidence = get_incidence(boundaries)
            % (neighbour, [element side], face) from the (side, element)
            % face table; the single neighbour of a boundary face is
            % duplicated
            [sorted_faces, p] = sort(boundaries(:)); % faces with multiplicity
            is_single = [diff(sorted_faces) ~= 0; true] & [true; diff(sorted_faces) ~= 0];
            p = repelem(p, 1 + is_single);           % boundary faces: duplicate incidence
            [s, K] = ind2sub(size(boundaries), p);
            incidence = [reshape(K, 2, 1, []), reshape(s, 2, 1, [])];
        end
    end

    methods (Access=protected)
        function out = parenReference(obj,indexOp)
            ErrorForUnsupportedIndexing(indexOp);
            out = submesh(obj, indexOp.Indices{1});
        end

        function obj = parenAssign(obj,indexOp,in)
            ErrorForUnsupportedIndexing(indexOp);
            ErrorForUnexpectedRightHandSide(in);
            return
        end
        function n = parenListLength(~,~,~)
            n = 1;
        end

        function obj = parenDelete(obj,indexOp)
            obj = submesh(obj, setdiff(1:size(obj.elements, 2), indexOp.Indices{1}));
        end
    end

    methods (Access=public)
        function out = cat(~, varargin) %#ok<STOUT>
            error("SimplicialGrid:cat", "Concatenation of meshes is not supported.");
        end

        function varargout = size(obj, varargin)
            sz = [size(obj.elements, 2), size(obj.faces, 2), size(obj.coordinates, 2)];
            varargout = cell(max(nargout, 1), 1);
            if isempty(varargin)
                if nargout <= 1
                    varargout{1} = sz(1);
                else
                    sz_cell = num2cell(sz);
                    [varargout{:}] = deal(sz_cell{:});
                end
            else
                if isequal(varargin{1}, "all")
                    varargout{1} = sz;
                elseif length(varargin{1}) > 1
                    varargout{1} = sz(varargin{1});
                else
                    sz_cell = cellfun(@(i) sz(i), varargin, 'UniformOutput', 0);
                    [varargout{:}] = deal(sz_cell{:});
                end
            end
        end
    end
end

function ErrorForUnsupportedIndexing(indexOp)
if ~isscalar(indexOp)
    error('Indexing after parentheses indexing is not supported.');
end
end
