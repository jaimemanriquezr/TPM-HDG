function [P, dP] = viscous_piola_spatial(G, G_dot, viscosity, rate)
% Kelvin–Voigt first Piola stress of the spatial rate of deformation
% (Newtonian Cauchy stress 2μ_v D), pointwise:
%   P = 2μ_v J D F⁻ᵀ = μ_v J (L + Lᵀ) F⁻ᵀ,   L = Ḟ F⁻¹,  D = sym L,
%   F = I + G,  J = det F,  Ḟ = Ġ,
% and its derivative along Ġ = rate G + const, with A = F⁻¹ and H the unit
% gradient b (dJ = J tr(A H), dA = −A H A, dL = (rate H − L H) A),
%   dP(:, a, b) = ∂P_a/∂G_b
%               = tr(A H) P + μ_v J ((dL + dLᵀ) Aᵀ − (L + Lᵀ) Aᵀ Hᵀ Aᵀ).
% Layouts as in viscous_piola.
n = size(G, 1);
F = G + [1 0 0 1];
J = F(:, 1).*F(:, 4) - F(:, 2).*F(:, 3);
A = [F(:, 4), -F(:, 2), -F(:, 3), F(:, 1)]./J;
At = tr2(A);
L = mm(G_dot, A);
DA = mm(L + tr2(L), At);
P = viscosity*J.*DA;
dP = zeros(n, 4, 4);
E4 = eye(4);
for b = 1:4
    H = repmat(E4(b, :), n, 1);
    dL = mm(rate*H - mm(L, H), A);
    dP(:, :, b) = sum(A.*tr2(H), 2).*P ...
        + viscosity*J.*(mm(dL + tr2(dL), At) - mm(DA, mm(tr2(H), At)));
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
