function export_interface_traction(stem, fluid_solution, solid_solution, ...
    fluid_mesh, solid_mesh, curve, own_interface, options)
%-------------------------------------------------------------------------
%Description: Writes the traction along the interface as plain CSV.
%
%             <stem>_true.csv   the stress extrapolated to the quadrature
%                               points of curve and contracted with its
%                               normal; comparable across interface
%                               discretizations and against the exact
%                               traction when one is supplied.
%             <stem>_own.csv    the same on own_interface, the interface
%                               the solve used, with its normal. For a
%                               polygonal interface this is the
%                               chord-normal traction actually imposed;
%                               when own_interface is curve it coincides
%                               with _true.
%
%             Columns: x, y, nx, ny, tS_x, tS_y, tF_x, tF_y [, tE_x, tE_y]
%             with n the solid-side normal of the solver's T1 convention
%             (n_S = +quadrature normals, n_F = -n_S), so a balanced
%             interface has tS + tF = T1 (zero on a physical wall).
%             tE is the exact solid traction sigma_S(x) n when
%             options.exact_sigma_S is given.
%
%             Points are listed patch by patch in partition order; sort by
%             the curve parameter when plotting.
%
%See also: transfer_traction_to_interface, export_fields_csv
%-------------------------------------------------------------------------
arguments
    stem            (1,1)   string
    fluid_solution  (1,1)   FESolution
    solid_solution  (1,1)   FESolution
    fluid_mesh              SimplicialGrid
    solid_mesh              SimplicialGrid
    curve                   AbstractCurve   % the exact interface
    own_interface           AbstractCurve   % the interface the solve used
    options.Order (1,1) double
    % 2x2 cell of @(x, y) handles for the exact solid stress, or empty.
    options.exact_sigma_S = []
end
k = options.Order;
write_one(stem + "_true.csv", curve, k, fluid_solution, solid_solution, ...
    fluid_mesh, solid_mesh, options.exact_sigma_S);
write_one(stem + "_own.csv", own_interface, k, fluid_solution, solid_solution, ...
    fluid_mesh, solid_mesh, options.exact_sigma_S);
end

%=========================================================================%
function write_one(file, iface, k, fs, ss, fluid_mesh, solid_mesh, exact_sigma_S)
[pF, pS] = iface.get_partitions(fluid_mesh, solid_mesh, ...
    quadrature_patch = k + 1, quadrature_path = k + 2);
nS = pS.quadrature.normals;
[xS, tS] = transfer_traction_to_interface(ss, solid_mesh, pS, Order = k, normals = nS);
[~, tF]  = transfer_traction_to_interface(fs, fluid_mesh, pF, Order = k, ...
    normals = -pF.quadrature.normals);
% Both partitions are induced by the same curve, so their interface points
% must coincide patch by patch; a mismatch would pair tractions at
% different points.
gap = max(abs(pF.quadrature.interface - pS.quadrature.interface), [], "all");
assert(gap <= 1e-10, "export_interface_traction:pointMismatch", ...
    "Fluid and solid interface points differ by %.3e.", gap);
X = reshape(xS, 2, []);  N = reshape(nS, 2, []);
M = [X; N; reshape(tS, 2, []); reshape(tF, 2, [])];
if ~isempty(exact_sigma_S)
    s = @(i, j) exact_sigma_S{i, j}(X(1, :), X(2, :)) + 0*X(1, :);
    tE = [s(1, 1).*N(1, :) + s(1, 2).*N(2, :); s(2, 1).*N(1, :) + s(2, 2).*N(2, :)];
    M = [M; tE];
end
writematrix(M.', file);
end
