function quadrature_points = ...
    get_induced_partition_quadrature(interface, mesh, induced_partition, options)
    %GET_INDUCED_PARTITION_QUADRATURE  Patch and transfer-path quadrature points.
    %
    %   q = get_induced_partition_quadrature(interface, mesh, ...
    %       induced_partition, quadrature_patch = 1, quadrature_path = 2)
    %
    %   evaluates, for the patches of one mesh's induced partition (see
    %   get_partitions), Gauss-Legendre points along each patch
    %   (quadrature_patch per patch) and along each transfer path joining
    %   the point on the mesh to the point on interface Sigma
    %   (quadrature_path per path). The option get_quadrature_points,
    %   @(v, N) points between the endpoints v(:,1,:,:) and v(:,2,:,:),
    %   replaces the Gauss-Legendre rule.
    %
    %   q.local      1 x nq x n, patch parameter on the mesh face.
    %   q.physical   2 x nq x n, the mesh points.
    %   q.interface  2 x nq x n, the Sigma points (interface.get_2d_points).
    %   q.path       2 x nq x n x np, the points of each transfer path.
    %   q.normals    2 x nq x n, interface.unit_normal at the Sigma points.
    %   with nq = quadrature_patch, np = quadrature_path, n = patches.
    arguments
        interface
        mesh SimplicialGrid
        induced_partition
        options.quadrature_patch = 1;
        options.quadrature_path = 2;
        options.get_quadrature_points function_handle = @(v, N) v(:, 1, :, :) + (lgwt(N, 0, 1).').*(v(:, 2, :, :) - v(:, 1, :, :))
    end
    get_quadrature_points = options.get_quadrature_points;
    n_quadpoints_patch = options.quadrature_patch;
    n_quadpoints_segment = options.quadrature_path;

    % dims: 1 x 2 x n_of_patches
    s_endpoints = induced_partition.vertices;
    partition_faces = induced_partition.faces;

    % dims: 1 x n_of_gauss x n_of_patches
    s_quad = get_quadrature_points(s_endpoints, n_quadpoints_patch);

    % dims: 2 x n_of_gauss x n_of_patches
    x_quad = get_1d_physical(mesh, s_quad, partition_faces);

    % Use canonical interface points from the parametrization so that
    % both meshes in a coupled problem evaluate exactly the same
    % interface point at each patch quadrature node
    s_interface_quad = get_quadrature_points(induced_partition.vertices_interface, n_quadpoints_patch);
    x_interface_quad = interface.get_2d_points(s_interface_quad);

    % dims: 2 x n_of_gauss x n_of_patches x n_of_gauss_CS
    x_tpath_endpoints = cat(2, permute(x_interface_quad, [1 4 3 2]), permute(x_quad, [1 4 3 2]));
    x_tpath_quad = get_quadrature_points(x_tpath_endpoints, n_quadpoints_segment);
    
    quadrature_points.local = s_quad;
    quadrature_points.physical = x_quad;
    quadrature_points.interface = x_interface_quad;
    quadrature_points.path = permute(x_tpath_quad, [1 4 3 2]);

    x_interface_quad_local = s_interface_quad;

    % each curve owns its normal convention (rotated tangent on the
    % parametrized lineage, grad phi on the level-set one)
    quadrature_points.normals = double(interface.unit_normal(x_interface_quad_local));
end