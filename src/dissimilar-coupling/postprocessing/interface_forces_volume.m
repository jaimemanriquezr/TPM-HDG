function [F, parts] = interface_forces_volume(solution, mesh, options)
%-------------------------------------------------------------------------
%Description: Force on an immersed body by the volume (Babuska-Miller /
%             virtual-work) formula, from the fluid's stress-form HDG
%             solution alone -- no boundary traction, no transfer:
%
%               F_i = - int_{Omega_f} [ sigma_h : grad(psi e_i)
%                                      + rho ((u_h . grad) u_h) . (psi e_i) ] dx
%
%             psi is a radial C^1 cutoff about `center`: 1 for
%             r <= r_inner, a cubic decay to 0 at r_outer, 0 beyond; the
%             outer boundary and the far field carry no weight. With the
%             steady momentum balance the integral equals the loop
%             integral of the numerical traction on the body up to the
%             discretization error, whatever the shape of the carved
%             boundary. Drag = F(1), lift = F(2).
%             Valid when the body lies in the disk r <= r_inner about
%             center (psi = 1 on it), the disk r < r_outer contains no
%             other boundary of the fluid mesh, and the fluid source term
%             vanishes for r < r_outer.
%
%             sigma_h is the full stress (sigma = 2 mu eps(u) - p I,
%             p = -tr(sigma)/2), so the pressure and viscous parts are
%             sigma_p = tr(sigma)/2 I and sigma_v = sigma - sigma_p. The
%             convective term (density > 0) uses grad u = eps(u) + gamma_h,
%             eps(u) = sigma_v / (2 mu), gamma_h = gamma [0 1; -1 0] the
%             skew part from the scalar rotation unknown.
%
%             Options: Order (quadrature degree Order + 3; default the
%             solution's), center, r_inner, r_outer (the cutoff), density
%             (> 0 adds the convective term) and viscosity (used by that
%             term only). parts holds the pressure, viscous and convective
%             contributions; F is their sum.
%-------------------------------------------------------------------------
arguments
    solution (1,1) FESolution      % fluid solution: fields (sigma, u, gamma)
    mesh           SimplicialGrid  % the fluid mesh
    options.Order (1,1) double = solution.Order
    options.center (2,1) double = [0; 0]
    options.r_inner (1,1) double {mustBePositive} = 0.07
    options.r_outer (1,1) double {mustBePositive} = 0.15
    options.density (1,1) double {mustBeNonnegative} = 0   % > 0: add the convective term
    options.viscosity (1,1) double {mustBePositive} = 1
end
assert(options.r_outer > options.r_inner, "interface_forces_volume:radii", ...
    "r_outer must exceed r_inner.");
k = options.Order;
[xq, w] = get_2d_reference_quadrature(k + 3);
w = reshape(w, [], 1)/sum(w);                           % int_K f = |K| sum_q w_q f(x_q)
n_el = size(mesh.elements, 2);
nq = size(xq, 2);

S = reshape(solution.Unknowns.field(1).evaluate(xq, mesh), nq, 2, 2, n_el);
X = mesh.get_2d_physical(xq);                           % 2 x nq x n_el
dx = squeeze(X(1, :, :)) - options.center(1);
dy = squeeze(X(2, :, :)) - options.center(2);
r = hypot(dx, dy);
s = (r - options.r_inner)/(options.r_outer - options.r_inner);
dpsi = -6*s.*(1 - s)/(options.r_outer - options.r_inner);   % d psi / d r of 1 - 3s^2 + 2s^3
dpsi(s <= 0 | s >= 1) = 0;
gx = dpsi .* dx ./ max(r, eps);
gy = dpsi .* dy ./ max(r, eps);

tr = (S(:, 1, 1, :) + S(:, 2, 2, :))/2;
Sp = tr .* reshape(eye(2), 1, 2, 2);
Sv = S - Sp;
vol = reshape(mesh.volume, 1, []);
integrate = @(A) -sum(vol .* sum(w .* A, 1));
contract = @(T) [integrate(squeeze(T(:, 1, 1, :)).*gx + squeeze(T(:, 1, 2, :)).*gy); ...
                 integrate(squeeze(T(:, 2, 1, :)).*gx + squeeze(T(:, 2, 2, :)).*gy)];
parts.pressure = contract(Sp);
parts.viscous = contract(Sv);
parts.convective = [0; 0];
if options.density > 0
    psi = 1 - 3*s.^2 + 2*s.^3;  psi(s <= 0) = 1;  psi(s >= 1) = 0;
    U = reshape(solution.Unknowns.field(2).evaluate(xq, mesh), nq, 2, n_el);
    gam = reshape(solution.Unknowns.field(3).evaluate(xq, mesh), nq, 1, 1, n_el);
    G = zeros(nq, 2, 2, n_el);  G(:, 1, 2, :) = gam;  G(:, 2, 1, :) = -gam;
    gradu = Sv/(2*options.viscosity) + G;               % grad u = eps(u) + gamma
    conv = @(i) squeeze(gradu(:, i, 1, :)).*squeeze(U(:, 1, :)) + squeeze(gradu(:, i, 2, :)).*squeeze(U(:, 2, :));
    parts.convective = options.density*[integrate(conv(1).*psi); integrate(conv(2).*psi)];
end
F = parts.pressure + parts.viscous + parts.convective;
end
