function [mesh_1, mesh_2] = fitted_split(interface, h, gap, options)
%-------------------------------------------------------------------------
%Description: Splits the unit square into two subdomains meeting at
%             `interface`, meshing each as its own polygon so that its
%             boundary follows the interface. A gap is opened explicitly,
%             by offsetting the two polygons by `gap`.
%
%             The two returned meshes are always ordered
%                 mesh_1  the outer / upper side
%                 mesh_2  the inner / lower side.
%
%             Topology is not an argument: it is read off the interface,
%             via interface.is_periodic().
%
%             This is the counterpart of SimplicialGrid.split, which
%             carves two subdomains out of a background grid instead.
%
%Usage:       [upper, lower] = fitted_split(SmoothCurveHorizontal(), h)
%             [outer, inner] = fitted_split(CircularCurve(), h, h^2)
%
%See also: SimplicialGrid.split, generate_mesh_list
%-------------------------------------------------------------------------
arguments
    interface ParametrizedCurve = SmoothCurveHorizontal();

    % Mesh size. A scalar applies to both sides; [h_1, h_2] sizes the
    % outer/upper and inner/lower meshes independently.
    h double = .2;

    % Gap. The two polygons are offset by this much so that Sigma sits
    % between the discrete interfaces; 0 leaves them coincident. The offset
    % is vertical for an open interface and normal for a closed one (see
    % open_polygons and closed_polygons).
    gap double = 0;

    options.boundary_markers double = [1, 1];

    % "unstructured" meshes each polygon with mesh2d; "structured" maps a
    % StructuredSimplicialGrid onto each side and is available for the open
    % topology only.
    options.type (1,1) string ...
        {mustBeMember(options.type, ["unstructured", "structured"])} = "unstructured"

    % Closed topology only. Number of polygon points around the whole
    % interface, as a function of h.
    options.interface_points function_handle = @(h) round(1/h) + 1
end
if isscalar(h); h = [h, h]; end

if interface.is_periodic()
    [mesh_1, mesh_2] = closed_polygons(interface, h, gap, options);
else
    [mesh_1, mesh_2] = open_polygons(interface, h, gap, options);
end
end

%========================== FITTED, OPEN TOPOLOGY ========================%
function [upper_mesh, lower_mesh] = open_polygons(I, h, d, options)
% The interface is a graph y = f(x) spanning the square, so the gap is
% opened by displacing it vertically, keeping both polygons graphs over
% the same x.
bdy_marker_1 = options.boundary_markers(1);
bdy_marker_2 = options.boundary_markers(2);

N1 = round(1/h(1)) + 1;
% Rounding the product, not the ratio, keeps both directions meaningful:
% h(2) = 2*h(1) halves the points, h(2) = h(1)/2 doubles them.
N2 = round(h(1) / h(2) * N1);

switch options.type
    case "unstructured"
        t1 = linspace(0, 1, N1);
        upper_rectangle = [[t1; I.y_coord(t1)+d], [1; 1], [0; 1]];
        upper_interface_conn = [[(1:N1-1); (2:N1)]; 3*ones(1, N1-1)];
        upper_conn = [upper_interface_conn, ...
            [[N1; N1+1]; bdy_marker_1], ...
            [[N1+1; N1+2]; bdy_marker_1], ...
            [[N1+2; 1]; bdy_marker_1]];
        upper_domain = {upper_rectangle; upper_conn};
        upper_mesh = SimplicialGrid(upper_domain{:}, ...
            element_diameter = h(1), ...
            mesh_generator = "mesh2d");

        t2 = linspace(1, 0, N2);
        lower_rectangle = [[0; 0], [1; 0], [t2; I.y_coord(t2)-d]];
        lower_interface_conn = [[(3:N2+1); (4:N2+2)]; 3*ones(1, N2-1)];
        lower_conn = [[[1; 2]; bdy_marker_2], ...
            [[2; 3]; bdy_marker_2], ...
            lower_interface_conn, ...
            [[N2+2; 1]; bdy_marker_2]];
        lower_domain = {lower_rectangle; lower_conn};
        lower_mesh = SimplicialGrid(lower_domain{:}, ...
            element_diameter = h(2), ...
            mesh_generator = "mesh2d");
    case "structured"
        upper_mesh = StructuredSimplicialGrid(N1, round(N1/2));
        y0 = I.y_coord(0);
        upper_mesh = [1, 0; 0 (y0-d)]*upper_mesh + [0; y0 + d];
        upper_mesh = apply_tag(upper_mesh, "Interface", @(X) isapprox(X(2, :), y0+d));

        lower_mesh = StructuredSimplicialGrid(N2, round(N2/2));
        lower_mesh = [1, 0; 0 (y0-d)]*lower_mesh;
        lower_mesh = apply_tag(lower_mesh, "Interface", @(X) isapprox(X(2, :), y0-d));
end
end

%========================= FITTED, CLOSED TOPOLOGY =======================%
function [outer_mesh, inner_mesh] = closed_polygons(I, h, d, options)
% A closed interface is not a graph over either axis, so here the gap is
% opened along the interface normal, leaving Sigma between the polygons.
% A valid offset only while d stays below the local radius of curvature.
if options.type ~= "unstructured"
    error("fitted_split:structuredClosed", ...
        "A structured fitted mesh is available for the open topology only.");
end
bdy_marker = options.boundary_markers(1);

N1 = options.interface_points(h(1));
t1 = linspace(I.t_range(1), I.t_range(2), N1+1);
t1 = t1(1:end-1);

vertices = [[0; 0], [1; 0], [1; 1], [0; 1]];
outer_rectangle = [vertices, [I.x_coord(t1); I.y_coord(t1)] + d*outward_normal(I, t1)];
outer_interface_conn = [4 + [(1:N1); (2:N1), 1]; 3*ones(1, N1)];
outer_conn = [[1; 2; bdy_marker], ...
    [2; 3; bdy_marker], ...
    [3; 4; bdy_marker], ...
    [4; 1; bdy_marker], ...
    outer_interface_conn];
outer_domain = {outer_rectangle; outer_conn};

outer_mesh = SimplicialGrid(outer_domain{:}, ...
    element_diameter = h(1), ...
    mesh_generator = "mesh2d");

N2 = options.interface_points(h(2));
t2 = linspace(I.t_range(1), I.t_range(2), N2+1);
t2 = t2(1:end-1);

inner_rectangle = [I.x_coord(t2); I.y_coord(t2)] - d*outward_normal(I, t2);
inner_interface_conn = [[(1:N2); (2:N2), 1]; 3*ones(1, N2)];
inner_domain = {inner_rectangle; inner_interface_conn};
inner_mesh = SimplicialGrid(inner_domain{:}, ...
    element_diameter = h(2), ...
    mesh_generator = "mesh2d", ...
    refine_boundaries = true);
end

function n = outward_normal(I, t)
% Tangent rotated by -90 degrees. Same convention as
% prune_infeasible_elements: for a counter-clockwise parametrization this
% points away from the enclosed region.
gamma_prime = I.tangent(t);
n = [gamma_prime(2, :); -gamma_prime(1, :)] ./ vecnorm(gamma_prime, 2, 1);
end
