function xhat = reference_element(dim, element_type)
%REFERENCE_ELEMENT  Vertices of the reference element.
%
%   xhat = reference_element(dim, element_type)
%
%   dim = 1 (segment) returns [0, 1]; dim = 2 (default) returns, as
%   columns, the vertices of the reference triangle (1, 0), (0, 1), (0, 0)
%   for element_type "Simplex" (default), or of the reference square
%   (0, 0), (1, 0), (1, 1), (0, 1) for "Quad", in the vertex order of the
%   elements of SimplicialGrid and CartesianGrid; any other dim returns [].
    arguments
        dim  double {mustBeScalarOrEmpty} = 2
        element_type (1,1) string {mustBeMember(element_type, ["Simplex", "Quad"])} = "Simplex"
    end
    switch dim
        case 1
            xhat = [0, 1];
        case 2
            switch element_type
                case "Simplex"
                    xhat = [1 0 0; 0 1 0];
                case "Quad"
                    xhat = [0 1 1 0; 0 0 1 1];
            end
        otherwise
            xhat = [];
    end
end