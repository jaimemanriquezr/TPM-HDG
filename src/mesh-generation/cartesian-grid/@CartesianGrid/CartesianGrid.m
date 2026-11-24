classdef CartesianGrid < AbstractGrid
    %CartesianGrid  Tensor-product grid of axis-aligned rectangles.
    %
    %   mesh = CartesianGrid(x, y)   breakpoints x(1) < ... < x(Nx+1) and
    %                                y(1) < ... < y(Ny+1), not necessarily
    %                                uniform
    %   mesh = CartesianGrid(N, M)   uniform N x M grid of the unit square
    %   mesh = CartesianGrid(N)      uniform N x N grid of the unit square
    %   mesh = CartesianGrid()       the empty mesh (CartesianGrid.empty)
    %
    %   Sibling of SimplicialGrid for Q^k discretizations (QkBasis); it is
    %   not accepted by FEProblem/solve. It stores the same tables as
    %   SimplicialGrid, with 4 sides per element, plus the breakpoints
    %   x, y. All numbering is lexicographic with the x index fastest, so
    %   data indexed by elements or faces reshape directly into matrices:
    %
    %     vertex (i, j)            i + (Nx+1)*(j-1),  i <= Nx+1, j <= Ny+1
    %     element (i, j)           i + Nx*(j-1),      i <= Nx,   j <= Ny
    %     vertical face (i, j)     i + (Nx+1)*(j-1),  i <= Nx+1, j <= Ny
    %       x = x(i), from y(j) to y(j+1)
    %     horizontal face (i, j)   n_v + i + Nx*(j-1), i <= Nx,  j <= Ny+1
    %       y = y(j), from x(i) to x(i+1); n_v = (Nx+1)*Ny
    %
    %   vertical_faces ((Nx+1) x Ny) and horizontal_faces (Nx x (Ny+1))
    %   hold these face indices as matrices.
    %
    %   Reference element [0, 1]^2, mapped by x = x_SW + diag(hx, hy) * r.
    %   Local sides: 1 south (s = 0), 2 east (r = 1), 3 north (s = 1),
    %   4 west (r = 0). Every face is parametrized in the +x (horizontal)
    %   or +y (vertical) direction, and so is every local side, so both
    %   elements of a face see the same parameter and permutation is all
    %   ones. As in SimplicialGrid, normals holds the face normal
    %   [dy; -dx] of each side's face (v1 -> v2, scaled by the length),
    %   +x for vertical and -y for horizontal faces, and the outward
    %   normal is orientation .* normals, so orientation = [1; 1; -1; -1].
    %
    %   Boundary faces are tagged "Dirichlet"; interior faces carry tag 0.
    %   Tag queries and apply_tag are inherited from AbstractGrid.
    %   submesh (and so mesh(element_list)) is not supported: a subset of
    %   the elements is no longer a tensor-product grid.
    %
    %See also: AbstractGrid, SimplicialGrid, QkBasis
    properties (Constant)
        geometric_dimension = 2
        ElementType = "Quad"
    end

    properties (Hidden)
        x               (1,:) double                    % breakpoints in x
        y               (1,:) double                    % breakpoints in y
    end

    properties
        coordinates     (2,:) double                    % physical coord.
        elements        (4,:) double                    % (vertex SW SE NE NW, element)
        faces           (3,:) double                    % (vertex 1, vertex 2, tag; face)
        boundaries      (4,:) double                    % (face S E N W, element)
        incidence       (2,2,:) double                  % (neighbour, [element side], face)
        orientation     (4,:) double                    % integer -1, 1
        permutation     (4,:) double                    % all ones
        normals         (8,:) double                    % face normals by side
        volume          (1,:) double                    % element areas
        area            (1,:) double                    % face lengths
        tau             (4,:) double                    % HDG stabilization
    end

    properties (Hidden, Dependent)
        Nx                                             % elements in x
        Ny                                              % elements in y
        hx              (1,:) double                    % diff(x)
        hy              (1,:) double                    % diff(y)
        vertical_faces                                  % (Nx+1) x Ny face indices
        horizontal_faces                                % Nx x (Ny+1) face indices
    end

    methods
        function obj = CartesianGrid(x, y)
            arguments
                x (1,:) double = []
                y (1,:) double = x
            end
            if isempty(x)
                return
            end
            if isscalar(x)
                x = linspace(0, 1, x + 1);
            end
            if isscalar(y)
                y = linspace(0, 1, y + 1);
            end
            assert(numel(x) >= 2 && all(diff(x) > 0), ...
                "CartesianGrid:breakpoints", "x must be strictly increasing.")
            assert(numel(y) >= 2 && all(diff(y) > 0), ...
                "CartesianGrid:breakpoints", "y must be strictly increasing.")
            obj.x = x;
            obj.y = y;
            Nx = numel(x) - 1;
            Ny = numel(y) - 1;

            [X, Y] = ndgrid(x, y);
            obj.coordinates = [X(:).'; Y(:).'];

            % elements: SW, SE, NE, NW (counter-clockwise)
            vtx = reshape(1:(Nx+1)*(Ny+1), Nx+1, Ny+1);
            sw = vtx(1:Nx, 1:Ny);
            obj.elements = [sw(:), sw(:) + 1, sw(:) + Nx + 2, sw(:) + Nx + 1].';

            obj.faces = CartesianGrid.get_faces(obj.elements);
            [obj.boundaries, obj.orientation, obj.incidence, obj.permutation] = ...
                CartesianGrid.get_boundaries(obj.faces, obj.elements);

            hx = diff(x);
            hy = diff(y);
            obj.volume = kron(hy, hx);
            obj.area = [kron(hy, ones(1, Nx+1)), repmat(hx, 1, Ny+1)];
            d = obj.coordinates(:, obj.faces(2, :)) - obj.coordinates(:, obj.faces(1, :));
            n_by_face = [d(2, :); -d(1, :)];
            obj.normals = reshape(n_by_face(:, obj.boundaries), 2*4, []);
            obj.tau = ones(4, Nx*Ny);
        end

        function n = get.Nx(obj)
            n = numel(obj.x) - 1;
        end

        function n = get.Ny(obj)
            n = numel(obj.y) - 1;
        end

        function h = get.hx(obj)
            h = diff(obj.x);
        end

        function h = get.hy(obj)
            h = diff(obj.y);
        end

        function F = get.vertical_faces(obj)
            F = reshape(1:(obj.Nx+1)*obj.Ny, obj.Nx+1, obj.Ny);
        end

        function F = get.horizontal_faces(obj)
            F = (obj.Nx+1)*obj.Ny + reshape(1:obj.Nx*(obj.Ny+1), obj.Nx, obj.Ny+1);
        end

        function new_obj = submesh(~, ~) %#ok<STOUT>
            error("CartesianGrid:submesh", ...
                "A subset of a CartesianGrid is not a tensor-product grid; submesh is not supported.");
        end
    end

    methods (Static)
        function out = empty()
            out = CartesianGrid();
        end

        function faces = get_faces(elements)
            % faces (vertex 1, vertex 2, tag) of the lexicographic grid
            % whose elements (SW, SE, NE, NW) are given: vertical faces
            % (+y), then horizontal faces (+x); boundary faces tagged 1
            % ("Dirichlet"), interior faces 0
            arguments
                elements (4, :)
            end
            [Nx, Ny] = grid_size(elements);
            vtx = reshape(1:(Nx+1)*(Ny+1), Nx+1, Ny+1);
            fv = vtx(:, 1:Ny);
            fh = vtx(1:Nx, :);
            on_bdry_v = false(Nx+1, Ny);  on_bdry_v([1, end], :) = true;
            on_bdry_h = false(Nx, Ny+1);  on_bdry_h(:, [1, end]) = true;
            faces = [fv(:).', fh(:).'; ...
                fv(:).' + Nx + 1, fh(:).' + 1; ...
                on_bdry_v(:).', on_bdry_h(:).'];
        end

        function [boundaries, orientation, incidence, permutation] = get_boundaries(faces, elements)
            % faces of the sides S, E, N, W of each element, their
            % orientation, the incidence table and the (all ones)
            % permutation, for faces numbered as in get_faces
            arguments
                faces (3, :)
                elements (4, :)
            end
            [Nx, Ny] = grid_size(elements);
            n_fv = (Nx+1)*Ny;
            assert(size(faces, 2) == n_fv + Nx*(Ny+1), ...
                "CartesianGrid:faces", "faces are not numbered as in get_faces.")
            Fv = reshape(1:n_fv, Nx+1, Ny);
            Fh = n_fv + reshape(1:Nx*(Ny+1), Nx, Ny+1);
            S = Fh(:, 1:Ny);  N = Fh(:, 2:Ny+1);
            W = Fv(1:Nx, :);  E = Fv(2:Nx+1, :);
            boundaries = [S(:), E(:), N(:), W(:)].';
            orientation = repmat([1; 1; -1; -1], 1, Nx*Ny);
            incidence = CartesianGrid.get_incidence(boundaries);
            permutation = ones(4, Nx*Ny);
        end
    end
end

function [Nx, Ny] = grid_size(elements)
% grid dimensions from the lexicographic connectivity: NW - SW = Nx + 1
Nx = elements(4, 1) - elements(1, 1) - 1;
Ny = size(elements, 2) / Nx;
end
