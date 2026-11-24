function [F, M] = interface_loop_resultants(t, x, partition, x0, line_element)
%-------------------------------------------------------------------------
%Description: Loop resultants of a traction field on the interface: net
%             force F = int t ds and torque M = int (x - x0) x t ds,
%             integrated over the induced partition's interface patches
%             with the patch quadrature the partition carries.
%-------------------------------------------------------------------------
arguments (Input)
    t         (2,:,:) double  % traction at the interface points, 2 x s x n_patches
    x         (2,:,:) double  % the interface points, same layout
    partition (1,1)   struct  % induced partition (get_partitions)
    x0        (2,1)   double  % torque center
    % Line element per patch (1 x n_patches). Empty (the default) uses
    % partition.area, the patch length on the mesh's own interface face;
    % pass interface_patch_arc(partition, curve) for the arc length on Sigma.
    line_element double = []
end
arguments (Output)
    F (2,1) double
    M (1,1) double
end
s_quad = size(partition.quadrature.local, 2);
[~, wS] = lgwt(s_quad, 0, 1);
ds = line_element;
if isempty(ds), ds = partition.area(1, :, 1); end
lw = reshape(ds, 1, 1, []) .* reshape(wS, 1, []);
F = [sum(lw .* t(1, :, :), "all"); sum(lw .* t(2, :, :), "all")];
r = x - x0;
M = sum(lw .* (r(1, :, :) .* t(2, :, :) - r(2, :, :) .* t(1, :, :)), "all");
end
