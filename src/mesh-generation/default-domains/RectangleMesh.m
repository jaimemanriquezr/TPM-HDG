classdef RectangleMesh < SimplicialGrid
    % Mesh of the rectangle [a, b] x [c, d] with element diameter h.
    %
    %   mesh = RectangleMesh(a, b, c, d, h, bottom_side = "Dirichlet", ...
    %       right_side = ..., top_side = ..., left_side = ...)
    %
    % Defaults are [-1, 1] x [-1, 1] and h = 0.25. Each side option is a
    % boundary tag name of the SimplicialGrid (an unknown name raises
    % RectangleMesh:unknownTag).
    %
    % See also: UnitSquareMesh, SimplicialGrid
    methods
        function obj = RectangleMesh(a, b, c, d, h, options)
            arguments
                a double = -1
                b double = 1
                c double = -1
                d double = 1
                h double = .25

                options.bottom_side = "Dirichlet";
                options.right_side = "Dirichlet";
                options.top_side = "Dirichlet";
                options.left_side = "Dirichlet";
            end
            bm_dict = SimplicialGrid().boundary_tags;
            sides = [options.bottom_side, options.right_side, options.top_side, options.left_side];
            unknown = sides(~isKey(bm_dict, sides));
            if ~isempty(unknown)
                error("RectangleMesh:unknownTag", "Unknown boundary tag ""%s"" (known: %s).", ...
                    unknown(1), strjoin(keys(bm_dict), ", "));
            end
            vertices = [a b b a;c c d d];
            bm = bm_dict(sides);
            connectivity = [1:4; (2:4), 1; bm];

            obj = obj@SimplicialGrid(vertices, connectivity, ...
                element_diameter = h)
        end
    end
end
