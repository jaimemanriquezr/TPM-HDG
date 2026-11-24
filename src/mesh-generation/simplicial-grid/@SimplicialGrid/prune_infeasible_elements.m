function elements = prune_infeasible_elements(bg_grid, elements, interface, max_deviation, options)
%-------------------------------------------------------------------------
%Description: Feasibility filter for an active mesh carved out of a
%             background grid. Drops every element contributing an
%             interface face whose outward normal deviates from the exact
%             outward normal, at the point the face maps to, by more than
%             max_deviation degrees.
%
%             At 90 degrees the face projects to zero interface measure:
%             it carries no transmission condition, so its trace degrees
%             of freedom are unconstrained and the coupled system is
%             singular. Smaller bounds additionally flatten the staircase,
%             turning a protruding step into a face roughly parallel to
%             the interface, at the cost of a larger gap.
%
%             The test is signed, so a face pointing back into its own
%             subdomain (a "fin") fails it. The projection is bracketed
%             against a dense sampling of Sigma rather than taken from the
%             curve's closest-point query, which can land on the opposite
%             surface of a thin geometry (e.g. an airfoil).
%
%             Removing an element exposes new boundary faces, which may
%             themselves be infeasible, so this iterates to a fixed point.
%-------------------------------------------------------------------------
arguments (Input)
    bg_grid                 SimplicialGrid          % background grid
    elements                double                  % candidate element list
    interface       (1,1)   AbstractCurve
    max_deviation   (1,1)   double {mustBeInRange(max_deviation, 0, 90)}

    % Only faces within this distance of Sigma, in multiples of the
    % face's own length, are judged. submesh tags every exposed face
    % "Interface"; without this bound, faces exposed by the pruning far
    % from Sigma are judged too and the erosion runs across the domain.
    options.interface_band (1,1) double {mustBePositive} = 2
end
arguments (Output)
    elements double  % feasible subset, in the original order
end
elements = elements(:);

% Reference sampling of Sigma for the bracketed projection, built once.
t_ref = linspace(interface.t_range(1), interface.t_range(2), 8001);
if is_periodic(interface)
    t_ref(end) = [];        % the seam point repeats the first
end
P_ref = interface.get_2d_points(t_ref);
dt_ref = t_ref(2) - t_ref(1);

% Each sweep removes at least one element or returns, so numel sweeps can
% exhaust the list and one more is needed to observe the fixed point (it
% also covers the empty input, which is already converged).
for sweep = 0:numel(elements)
    mesh = submesh(bg_grid, elements.');
    faces = mesh.intfaces;
    if isempty(faces), return; end

    P1 = mesh.coordinates(1:2, mesh.faces(1, faces));
    P2 = mesh.coordinates(1:2, mesh.faces(2, faces));
    mid = (P1 + P2) / 2;

    edge = P2 - P1;
    normal_h = [edge(2, :); -edge(1, :)] ./ vecnorm(edge, 2, 1);
    % Orient each normal out of the element that owns the face: the
    % vertex order in mesh.faces is arbitrary, the owning centroid is not.
    [~, slot] = ismember(faces, mesh.boundaries(:));
    owner = ceil(slot / 3);
    centroid = squeeze(mean(reshape( ...
        mesh.coordinates(:, mesh.elements(:, owner)), 2, 3, []), 2));
    inward = sum(normal_h .* (mid - centroid), 1) < 0;
    normal_h(:, inward) = -normal_h(:, inward);

    t_mid = bracketed_parameter(interface, mid, t_ref, P_ref, dt_ref);
    normal_exact = interface.unit_normal(t_mid);

    % Only faces on the interface are judged -- see interface_band.
    on_interface = vecnorm(mid - interface.get_2d_points(t_mid)) ...
        <= options.interface_band * mesh.area(faces);
    if ~any(on_interface), return; end

    % unit_normal points out of the enclosed region, and this runs for
    % both subdomains, so the expected sign is read off the faces
    % themselves rather than passed in -- from the interface ones, whose
    % orientation is the only one that means anything.
    dotted = sum(normal_h .* normal_exact, 1);
    side = sign(median(dotted(on_interface)));
    if side == 0, side = 1; end

    alignment = side * dotted;
    infeasible = faces(on_interface & (alignment < cosd(max_deviation)));
    if isempty(infeasible), return; end

    % submesh preserves element order, so a local index addresses the same
    % entry of `elements`.
    local_elements = mesh.get_element_by_face(infeasible);
    elements(unique(local_elements)) = [];
end
error("prune_infeasible_elements:didNotConverge", ...
    "Feasibility pruning did not reach a fixed point.");
end

%=========================================================================%
function t = bracketed_parameter(interface, X, t_ref, P_ref, dt_ref)
% Parameter of the closest point of Sigma: nearest reference sample, then
% a ternary search on |X - gamma(t)|^2 inside that bracket. Derivative-free
% and vectorised: it needs only get_2d_points.
n = size(X, 2);
nearest = zeros(1, n);
chunk = 200;                        % bounds the distance block's memory
for a = 1:chunk:n
    b = min(a + chunk - 1, n);
    d2 = (P_ref(1, :).' - X(1, a:b)).^2 + (P_ref(2, :).' - X(2, a:b)).^2;
    [~, nearest(a:b)] = min(d2, [], 1);
end
lo = t_ref(nearest) - dt_ref;
hi = t_ref(nearest) + dt_ref;
squared_distance = @(t) sum((X - interface.get_2d_points(t)).^2, 1);
for iteration = 1:20               % bracket shrinks by (2/3)^20
    m1 = lo + (hi - lo) / 3;
    m2 = hi - (hi - lo) / 3;
    closer = squared_distance(m1) < squared_distance(m2);
    hi(closer)   = m2(closer);
    lo(~closer)  = m1(~closer);
end
t = (lo + hi) / 2;
end
