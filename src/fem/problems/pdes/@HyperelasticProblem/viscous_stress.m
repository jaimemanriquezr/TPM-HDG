function [P, dP] = viscous_stress(G, G_dot, viscosity, rate, model)
% Kelvin–Voigt first Piola stress and its derivative along
% Ġ = rate G + const, pointwise, F = I + G, Ḟ = Ġ, μ_v = viscosity:
%   model = "material"  P = F S_v,  S_v = 2μ_v Ė = μ_v(ḞᵀF + FᵀḞ)
%                       (viscous_piola)
%   model = "spatial"   P = 2μ_v J D F⁻ᵀ,  D = sym(Ḟ F⁻¹)
%                       (viscous_piola_spatial)
% G, G_dot are n x 4 (column-major [G11 G21 G12 G22], as in first_piola);
% P is n x 4 and dP(:, a, b) = ∂P_a/∂G_b, n x 4 x 4, includes the rate.
switch model
    case "material"
        [P, dP] = viscous_piola(G, G_dot, viscosity, rate);
    case "spatial"
        [P, dP] = viscous_piola_spatial(G, G_dot, viscosity, rate);
    otherwise
        error("HyperelasticProblem:viscousModel", ...
            'viscous_model must be "material" or "spatial", not "%s".', model);
end
end
