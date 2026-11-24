function [offset, e_rest] = fluid_stress_offset(fs, sigma_exact)
%FLUID_STRESS_OFFSET  Constant pressure part of the fluid stress error.
%
%   [offset, e_rest] = fluid_stress_offset(fs, sigma_exact)
%
%   fs: a fluid FESolution (field 1 the stress), sigma_exact: the 2 x 2
%   cell of the exact stress. With the pressure level fixed by the
%   coupling (pressure_mean = "pinned"), a global constant in the pressure
%   error shows as offset*I in sigma_h - sigma: offset is the mean over the
%   fluid mesh of tr(sigma_h - sigma)/2, e_rest the L2 norm of
%   sigma_h - sigma - offset*I (the error modulo that constant).
mesh = fs.Mesh;
field = fs.Unknowns.field(1);
basis = field.basis;
[xy, rs, w] = quadpoints_elements_2d(mesh, basis.Order + 1, output="array");
values = pagemtimes(basis.evaluate(rs, mesh), reshape(field.dof, basis.Dimension, 1, 1, 1, []));
n_q = size(values, 1);
E = cell(2, 2);
for i = 1:4
    [r, c] = ind2sub([2 2], i);
    sigma_h = reshape(values(:, 1, r, c, :), n_q, []);
    E{r, c} = sigma_h - reshape(sigma_exact{r, c}(xy(1, :, :), xy(2, :, :)) + 0*xy(1, :, :), size(sigma_h));
end
integral = @(X) sum(mesh.volume(:).' .* (w.'*X));
offset = integral((E{1, 1} + E{2, 2})/2)/sum(mesh.volume);
E{1, 1} = E{1, 1} - offset;
E{2, 2} = E{2, 2} - offset;
e_rest = sqrt(sum(cellfun(@(X) integral(X.^2), E), "all"));
end
