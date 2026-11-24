classdef StructuredSimplicialGrid < SimplicialGrid
    %STRUCTUREDSIMPLICIALGRID  Structured triangulation of a rectangle.
    %
    %   mesh = StructuredSimplicialGrid(N, M, a, b, c, d, style = ...)
    %   splits the rectangle (a, b) x (c, d) into N x M cells, each cut
    %   into triangles. Defaults: N = M = 10, (a, b) = (0, 1),
    %   (c, d) = (0, 1), style = "left". A non-integer N is read as a mesh
    %   size h and gives round(1/N) cells in each direction.
    %
    %   mesh = StructuredSimplicialGrid(N, M, style) is the short form on
    %   the unit square.
    %
    %   style selects how each cell is cut:
    %     "left"     two triangles, sharing the diagonal from the south-east
    %                to the north-west corner
    %     "right"    two triangles, sharing the diagonal from the south-west
    %                to the north-east corner
    %     "crossed"  four triangles around an added vertex at the cell
    %                centre
    %
    %   Every boundary face carries the tag "Dirichlet"; interior faces are
    %   untagged.
    %
    %See also: SimplicialGrid
    methods
        function obj = StructuredSimplicialGrid(N, M, a, b, c, d, options)
            arguments
                N = 10
                M = 10

                a = 0;
                b = 1;
                c = 0;
                d = 1;

                options.style = "left";
            end
            if isstring(a) || ischar(a)
                % (N, M, style) call form
                options.style = string(a);
                a = 0;
            end
            if rem(N, 1) > 0
                [N, M] = deal(round(1 / N));
            end

            x = linspace(a, b, N+1);
            y = linspace(c, d, M+1);

            [x_grid, y_grid] = meshgrid(x, y);
            input_vertices = [x_grid(:), y_grid(:)].';

            SW = (1:M) + permute((0:N-1)*(M+1), [1 3 2]);
            SE = (1:M) + (M+1) + permute((0:N-1)*(M+1), [1 3 2]);
            NE = (2:M+1) + (M+1) + permute((0:N-1)*(M+1), [1 3 2]);
            NW = (2:M+1) + permute((0:N-1)*(M+1), [1 3 2]);
            switch options.style
                case "left"
                    all_simplex = {
                        [SE; NW; SW];
                        [NW; SE; NE]
                        };
                case "right"
                    all_simplex = {
                        [NE; NW; SW];
                        [NE; SW; SE]
                        };
                case "crossed"
                    x_mid = mean(x([1:end-1;2:end]));
                    y_mid = mean(y([1:end-1;2:end]));
                    [x_grid_mid, y_grid_mid] = meshgrid(x_mid, y_mid);
                    mid_vertices = [x_grid_mid(:), y_grid_mid(:)].';
                    input_vertices = [input_vertices, mid_vertices];

                    C = (N+1)*(M+1) + (1:M) + permute((0:N-1)*M, [1 3 2]);

                    all_simplex = {[NW; C; NE];
                        [SW; C; NW];
                        [SE; C; SW];
                        [NE; C; SE]};
            end
            input_connectivity = reshape(cat(2, all_simplex{:}), 3, []);

            obj  = obj@SimplicialGrid(input_vertices, input_connectivity,...
                mesh_generator="array");
        end
    end
end
