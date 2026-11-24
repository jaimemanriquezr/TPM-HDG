function [upper_mesh, lower_mesh, I] = blood_vessel(interface, h, d, options)
%BLOOD_VESSEL  Pair of meshes of a vessel wall (upper) and lumen (lower) on [0, 6].
%
%   [upper_mesh, lower_mesh, I] = blood_vessel(interface, h, d, ...
%       interface_points = [])
%
%   The upper mesh lies between the interface (shifted up by d) and the
%   line y = 0.6, the lower one between y = 0 and the interface (shifted
%   down by d), both over x in [0, 6]; the interface edges are tagged 3
%   and the remaining edges 1. The default interface is the flat wall
%   y = 0.5 (a FlatCurveHorizontal with t_range [0 6]); I is the
%   interface curve used. interface, h, d and interface_points are
%   described at the arguments below.
    arguments
        % The interface, in the split_fun slot of generate_mesh_list.
        % Empty gives the flat wall y = 0.5; a curve (e.g. a stenosed
        % wall) must be a graph over [0, 6], since the gap d is applied
        % vertically and the solid's outer boundary stays at y = 0.6. A
        % constriction of depth delta narrows the lumen and thickens the
        % wall to 0.1 + delta, so delta is not bounded by the wall
        % thickness.
        interface = []
        % Mesh size. A scalar applies to both sides; [h_1, h_2] sizes the
        % upper and lower meshes independently.
        h double = .25/2;
        % Retreat of each mesh from the interface: a scalar retreats both
        % meshes by d (upper at 0.5 + d, lower at 0.5 - d); a 1x2 vector
        % [d_upper d_lower] places them independently, e.g. [0 2*d] keeps
        % the upper (solid) mesh conforming with the same total gap 2*d.
        d (1,:) double {mustBeNonnegative, mustBeNonempty} = 0;

        % Points sampling the curve on the upper mesh. Empty uses
        % round(1/h) + 1, a chord every 6h over [0, 6], so a curved wall
        % stays a 6h polygon; pass round(6/h) + 1 to sample it at h.
        options.interface_points double {mustBeScalarOrEmpty} = []
    end
    if isscalar(h); h = [h, h]; end
    d_upper = d(1);  d_lower = d(end);
    if isempty(interface)
        I = FlatCurveHorizontal(0.5);
        I.t_range = [0 6];
    else
        I = interface;
    end

    N1 = round(1/h(1)) + 1;
    if ~isempty(options.interface_points); N1 = options.interface_points; end
    % Same form as fitted_split: round the product, not the ratio, so both
    % refinement directions stay meaningful.
    N2 = round(h(1) / h(2) * N1);

    t1 = linspace(0, 6, N1);

    upper_rectangle = [[t1; I.y_coord(t1)+d_upper], [6; 0.6], [0; 0.6]];

    upper_interface_conn = [[(1:N1-1); (2:N1)]; 3*ones(1, N1-1)];
    upper_conn = [upper_interface_conn, ...
                  [[N1; N1+1]; 1], ...
                  [[N1+1; N1+2]; 1], ...
                  [[N1+2; 1]; 1]];
    upper_domain = {upper_rectangle; upper_conn};
    upper_mesh = SimplicialGrid(upper_domain{:}, ...
            element_diameter = h(1), ...
            mesh_generator = "mesh2d");

    t2 = linspace(6, 0, N2);
    lower_rectangle = [[0; 0], [6; 0], [t2; I.y_coord(t2)-d_lower]];
    lower_interface_conn = [[(3:N2+1); (4:N2+2)]; 3*ones(1, N2-1)];
    lower_conn = [[[1; 2]; 1], ...
                  [[2; 3]; 1], ...
                  lower_interface_conn, ...
                  [[N2+2; 1]; 1]];
    lower_domain = {lower_rectangle; lower_conn};
    lower_mesh = SimplicialGrid(lower_domain{:}, ...
            element_diameter = h(2), ...
            mesh_generator = "mesh2d");
end
