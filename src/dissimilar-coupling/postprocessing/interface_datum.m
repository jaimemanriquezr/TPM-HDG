function [datum, mapping] = interface_datum(interface, fluid_mesh, partition_solid, coeffs, options)
%-------------------------------------------------------------------------
%Description: A field on Sigma given per patch by Legendre coefficients in
%             the solid face-local coordinate, as a datum handle on Sigma,
%             with the map of the fluid interface faces onto Sigma that
%             goes with it.
%
%               [datum, mapping] = interface_datum(interface, fluid_mesh, ...
%                   partition_solid, coeffs)
%
%             coeffs (2*d1 x n_patches, x-block over y-block) has the
%             layout of the T2 history h0 of InterfaceProblem.solve and of
%             trace_extension_coefficients: column e lives on patch e of
%             partition_solid (the solid side of get_partitions), in the
%             face-local coordinate s of that patch's solid face.
%
%             datum = {@(x, y) ...; @(x, y) ...} evaluates the field at
%             points x_bar of Sigma: the curve parameter t of x_bar
%             (map_points_to_interface) selects the patch whose parameter
%             interval contains it, t is mapped affinely to s on that
%             patch, and the Legendre polynomial is evaluated there.
%
%             mapping(X) (2 x n -> 2 x n) maps points X on the fluid
%             faces tagged "Interface" onto Sigma as the TPM does
%             (map_mesh_to_interface): the face-local coordinate of X is
%             mapped affinely to the curve parameter of the face's
%             segment. On a curved Sigma this differs inside a face from
%             the closest point by O(h^2); at the vertices and on a flat
%             Sigma they coincide. options.mapping is the fluid mesh's
%             vertex rule, as passed to get_partitions as mapping_1 (""
%             keeps the curve's own). Pass both to
%             TransferEssentialBoundaryConditions to impose the field as
%             an essential datum on the fluid interface faces.
%
%             tip_mapping must be the one the partition was built with
%             (get_partitions). Under "collapse" a collapsed fluid face
%             maps every point to its collapse point gamma(t*), and the
%             datum is read from the regular patches only (a collapsed
%             patch is the single point t*, where the regular patches
%             already give the value).
%
%See also: FluidStructureCoupling.initial_state, trace_extension_coefficients,
%          AbstractCurve.get_partitions, AbstractCurve.map_mesh_to_interface
%-------------------------------------------------------------------------
arguments
    interface           AbstractCurve
    fluid_mesh          SimplicialGrid
    partition_solid     (1,1) struct
    coeffs              (:,:) double
    options.tip_mapping (1,1) string ...
        {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
    options.mapping     (1,1) string = ""
end
d1 = size(coeffs, 1)/2;
k = d1 - 1;
t_patch = reshape(partition_solid.vertices_interface, 2, []);   % [t_lo; t_hi]
s_patch = reshape(partition_solid.vertices, 2, []);             % s at t_lo, t_hi
assert(size(coeffs, 2) == size(t_patch, 2), "interface_datum:size", ...
    "coeffs must have one column per patch (%d), got %d.", size(t_patch, 2), size(coeffs, 2));
% patches with an interval of Sigma (all of them unless faces collapsed)
regular = 1:size(t_patch, 2);
if isfield(partition_solid, "collapsed")
    regular = find(~partition_solid.collapsed);
end
% the datum is evaluated at points of Sigma, where every map returns the
% point's own parameter; the "normal" rule needs per-point normals there,
% so the closest point stands in for it
on_sigma = interface;
if isprop(on_sigma, "mapping") && on_sigma.mapping == "normal"
    on_sigma.mapping = "closest";
end
datum = {@(x, y) evaluate(x, y, 1); @(x, y) evaluate(x, y, 2)};

segments = interface.map_mesh_to_interface(fluid_mesh, ...
    tip_mapping = options.tip_mapping, mapping = options.mapping);
faces = fluid_mesh.intfaces(:).';
V1 = fluid_mesh.coordinates(1:2, fluid_mesh.faces(1, faces));
V2 = fluid_mesh.coordinates(1:2, fluid_mesh.faces(2, faces));
mapping = @(X) to_interface(X);

    function value = evaluate(x, y, component)
        t = reshape(on_sigma.map_points_to_interface([x(:).'; y(:).']), [], 1);
        % the patch whose parameter interval contains t (the nearest one at
        % the ends)
        e = regular(discretize(t, [-Inf, t_patch(1, regular(2:end)), Inf]));
        s = s_patch(1, e).' + (t - t_patch(1, e).')./(t_patch(2, e) - t_patch(1, e)).' ...
            .* (s_patch(2, e) - s_patch(1, e)).';
        rows = (component - 1)*d1 + (1:d1);
        value = sum(legendre_basis(k, s.') .* coeffs(rows, e).', 2);
        value = reshape(value, size(x)) + 0*x;
    end

    function X_bar = to_interface(X)
        X = reshape(X, 2, []);
        n = size(X, 2);
        % closest fluid interface face and the face-local coordinate there
        D = V2 - V1;
        L2 = sum(D.^2, 1);
        s = ((X(1, :).' - V1(1, :)).*D(1, :) + (X(2, :).' - V1(2, :)).*D(2, :))./L2;
        s = min(max(s, 0), 1);
        dist = hypot(X(1, :).' - V1(1, :) - s.*D(1, :), X(2, :).' - V1(2, :) - s.*D(2, :));
        [d_min, j] = min(dist, [], 2);
        assert(all(d_min <= 1e-8*sqrt(L2(j)).'), "interface_datum:offFace", ...
            "mapping expects points on the fluid interface faces (distance %g).", max(d_min));
        s = s(sub2ind([n, numel(faces)], (1:n).', j));
        % the segment of that face whose s-range holds s (two on a face
        % straddling the seam of a periodic curve)
        t = zeros(n, 1);
        for i = 1:n
            collapsed = find(segments.collapsed_face == faces(j(i)), 1);
            if ~isempty(collapsed)
                t(i) = segments.collapsed_t(collapsed);
                continue
            end
            q = find(segments.face == faces(j(i)));
            if numel(q) > 1
                lo = min(segments.s(:, q), [], 1);  hi = max(segments.s(:, q), [], 1);
                q = q(find(s(i) >= lo & s(i) <= hi, 1));
            end
            t(i) = segments.t(1, q) + (s(i) - segments.s(1, q)) ...
                / (segments.s(2, q) - segments.s(1, q)) * (segments.t(2, q) - segments.t(1, q));
        end
        X_bar = reshape(interface.get_2d_points(reshape(t, 1, [])), 2, []);
    end
end
