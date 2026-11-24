function [x_hat, w, x] = get_quadrature(mesh, k, options)
%GET_QUADRATURE  Reference quadrature points of a mesh and their images.
%
%   [x_hat, w, x] = get_quadrature(mesh, k, coordinates = "local")
%
%   For a mesh of geometric dimension 2, x_hat and w are the rule of
%   get_2d_reference_quadrature(k, mesh.ElementType) and x its image under
%   mesh.get_2d_physical. For dimension 1 nothing is returned.
    arguments
        mesh AbstractGrid
        k double = 0
        options.coordinates = "local"
    end
    
    switch mesh.geometric_dimension
        case 1
            return;
        case 2
            switch options.coordinates
                case "local"
                    [x_hat, w] = get_2d_reference_quadrature(k, mesh.ElementType);
                    x = get_2d_physical(mesh, x_hat);
            end
    end
end
        