function [P, dP] = first_piola(G, lambda, mu, material)
% First Piola–Kirchhoff stress of the displacement gradient G and its
% derivative, pointwise (plane strain).
%
%   [P, dP] = HyperelasticProblem.first_piola(G, lambda, mu, material)
%
%   G   n x 4, the 2x2 gradients stored column-major [G11 G21 G12 G22]
%       (the component order of the PEERS stress blocks)
%   P   n x 4, P(I + G) in the same layout
%   dP  n x 4 x 4, dP(:, a, b) = ∂P_a/∂G_b
% material: "saint_venant_kirchhoff" (P = F S, S = λ tr E I + 2μ E,
% E = (FᵀF − I)/2) or "neo_hookean" (P = μ(F − F⁻ᵀ) + λ log(J) F⁻ᵀ, which
% raises HyperelasticProblem:inverted when J ≤ 0), or
% "neo_hookean_eulerian": the Cauchy stress of the same law as a function
% of the Eulerian displacement gradient G = ∇ₓu (X = x − u, F = (I − G)⁻¹),
% σ = J⁻¹ P Fᵀ = j(μ(B − I) − λ log(j) I), j = det(I − G), B = F Fᵀ
% (raises HyperelasticProblem:inverted when j ≤ 0); P and dP are then σ
% and ∂σ/∂G. "saint_venant_kirchhoff_eulerian": likewise the Cauchy stress
% of the Saint Venant–Kirchhoff law, σ = j F S Fᵀ, F = (I − G)⁻¹.
arguments
    G (:, 4) double
    lambda (1,1) double
    mu (1,1) double
    material (1,1) string = "saint_venant_kirchhoff"
end
n = size(G, 1);
F = G + [1 0 0 1];
E4 = eye(4);
dP = zeros(n, 4, 4);
switch material
    case "saint_venant_kirchhoff"
        S = svk_stress(F, lambda, mu);
        P = mm(F, S);
        for b = 1:4
            H = repmat(E4(b, :), n, 1);
            dE = sym2(mm(tr2(F), H));
            dS = lambda*(dE(:, 1) + dE(:, 4)).*[1 0 0 1] + 2*mu*dE;
            dP(:, :, b) = mm(H, S) + mm(F, dS);
        end
    case "neo_hookean"
        J = F(:, 1).*F(:, 4) - F(:, 2).*F(:, 3);
        if any(J <= 0)
            error("HyperelasticProblem:inverted", ...
                "The deformation gradient has det F <= 0 at %d point(s).", nnz(J <= 0));
        end
        FiT = [F(:, 4), -F(:, 3), -F(:, 2), F(:, 1)]./J;   % F⁻ᵀ
        logJ = log(J);
        P = mu*(F - FiT) + lambda*logJ.*FiT;
        for b = 1:4
            H = repmat(E4(b, :), n, 1);
            FiT_HT_FiT = mm(mm(FiT, tr2(H)), FiT);
            dP(:, :, b) = mu*H + (mu - lambda*logJ).*FiT_HT_FiT ...
                + lambda*sum(FiT.*H, 2).*FiT;
        end
    case "neo_hookean_eulerian"
        a = [1 0 0 1] - G;                                 % F⁻¹
        j = a(:, 1).*a(:, 4) - a(:, 2).*a(:, 3);           % 1/J
        if any(j <= 0)
            error("HyperelasticProblem:inverted", ...
                "The displacement gradient has det(I - G) <= 0 at %d point(s).", nnz(j <= 0));
        end
        F = [a(:, 4), -a(:, 2), -a(:, 3), a(:, 1)]./j;
        B = mm(F, tr2(F));
        P = j.*(mu*(B - [1 0 0 1]) - lambda*log(j).*[1 0 0 1]);
        for b = 1:4
            H = repmat(E4(b, :), n, 1);
            trFH = sum(F.*tr2(H), 2);
            FHB = mm(mm(F, H), B);
            dP(:, :, b) = -trFH.*P + mu*j.*(FHB + tr2(FHB)) + lambda*j.*trFH.*[1 0 0 1];
        end
    case "saint_venant_kirchhoff_eulerian"
        a = [1 0 0 1] - G;                                 % F⁻¹
        j = a(:, 1).*a(:, 4) - a(:, 2).*a(:, 3);           % 1/J
        if any(j <= 0)
            error("HyperelasticProblem:inverted", ...
                "The displacement gradient has det(I - G) <= 0 at %d point(s).", nnz(j <= 0));
        end
        F = [a(:, 4), -a(:, 2), -a(:, 3), a(:, 1)]./j;
        S = svk_stress(F, lambda, mu);
        P = j.*mm(mm(F, S), tr2(F));
        for b = 1:4
            H = repmat(E4(b, :), n, 1);
            trFH = sum(F.*tr2(H), 2);                      % −dj/j
            dF = mm(mm(F, H), F);
            dE = sym2(mm(tr2(F), dF));
            dS = lambda*(dE(:, 1) + dE(:, 4)).*[1 0 0 1] + 2*mu*dE;
            X = mm(mm(dF, S), tr2(F));
            dP(:, :, b) = -trFH.*P + j.*(X + tr2(X) + mm(mm(F, dS), tr2(F)));
        end
    otherwise
        error("HyperelasticProblem:material", "Unknown material ""%s"".", material);
end
end

function S = svk_stress(F, lambda, mu)
% S = λ tr E I + 2μ E, E = (FᵀF − I)/2
E = (mm(tr2(F), F) - [1 0 0 1])/2;
S = lambda*(E(:, 1) + E(:, 4)).*[1 0 0 1] + 2*mu*E;
end

function C = mm(A, B)
% 2x2 products, column-major n x 4 layout
C = [A(:, 1).*B(:, 1) + A(:, 3).*B(:, 2), A(:, 2).*B(:, 1) + A(:, 4).*B(:, 2), ...
     A(:, 1).*B(:, 3) + A(:, 3).*B(:, 4), A(:, 2).*B(:, 3) + A(:, 4).*B(:, 4)];
end

function A = tr2(A)
A = A(:, [1 3 2 4]);
end

function A = sym2(A)
A = (A + tr2(A))/2;
end
