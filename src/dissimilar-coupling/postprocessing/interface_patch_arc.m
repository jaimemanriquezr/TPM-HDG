function arc = interface_patch_arc(partition, curve, options)
%-------------------------------------------------------------------------
%Description: Arc length on Sigma of every patch of an induced partition
%             (get_partitions), 1 x n_patches, from the patch endpoints in
%             the interface parameter (partition.vertices_interface).
%
%             partition.area is the patch length on the mesh's own
%             interface face (the carved polygon), which is what the
%             face-tested transmission rows need; a loop integral on
%             Sigma (interface_loop_resultants) needs this arc length
%             instead. Exact for a CircularCurve (R |Delta theta|); a
%             polyline of curve.get_2d_points otherwise.
%-------------------------------------------------------------------------
arguments
    partition (1,1) struct
    curve     (1,1) AbstractCurve
    options.n_samples (1,1) double {mustBeInteger, mustBePositive} = 201
end
t_lo = reshape(partition.vertices_interface(1, 1, :), 1, []);
t_hi = reshape(partition.vertices_interface(1, 2, :), 1, []);
if isa(curve, "CircularCurve")
    % x(t) = c + R (cos 2 pi t, sin 2 pi t)
    arc = 2*pi*curve.radius*abs(t_hi - t_lo);
    return;
end
n_p = numel(t_lo);
arc = zeros(1, n_p);
for e = 1:n_p
    X = curve.get_2d_points(linspace(t_lo(e), t_hi(e), options.n_samples));
    arc(e) = sum(vecnorm(diff(X, 1, 2), 2, 1));
end
end
