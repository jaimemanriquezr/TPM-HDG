function [P, dP] = viscous_piola(G, G_dot, viscosity, rate)
% Kelvin–Voigt first Piola stress of the Green–Lagrange strain rate,
% pointwise:
%   P = F S,   S = 2μ_v Ė = μ_v (ḞᵀF + FᵀḞ),   F = I + G,  Ḟ = Ġ,
% and its derivative along Ġ = rate G + const,
%   dP(:, a, b) = ∂P_a/∂G_b = H S + F μ_v(rate HᵀF + ḞᵀH + HᵀḞ + rate FᵀH),
% H the unit gradient b. G, G_dot are n x 4 (column-major
% [G11 G21 G12 G22], as in first_piola); P is n x 4, dP n x 4 x 4.
n = size(G, 1);
F = G + [1 0 0 1];
S = viscosity*(mm(tr2(G_dot), F) + mm(tr2(F), G_dot));
P = mm(F, S);
dP = zeros(n, 4, 4);
E4 = eye(4);
for b = 1:4
    H = repmat(E4(b, :), n, 1);
    dS = viscosity*(rate*mm(tr2(H), F) + mm(tr2(G_dot), H) + mm(tr2(H), G_dot) ...
        + rate*mm(tr2(F), H));
    dP(:, :, b) = mm(H, S) + mm(F, dS);
end
end

function C = mm(A, B)
% 2x2 products, column-major n x 4 layout
C = [A(:, 1).*B(:, 1) + A(:, 3).*B(:, 2), A(:, 2).*B(:, 1) + A(:, 4).*B(:, 2), ...
     A(:, 1).*B(:, 3) + A(:, 3).*B(:, 4), A(:, 2).*B(:, 3) + A(:, 4).*B(:, 4)];
end

function A = tr2(A)
A = A(:, [1 3 2 4]);
end
