function P = first_piola_symbolic(G, lambda, mu, material)
% Symbolic first Piola–Kirchhoff stress P(I + G) of a 2x2 symbolic
% displacement gradient G (the laws of HyperelasticProblem.first_piola).
arguments
    G (2, 2) sym
    lambda
    mu
    material (1,1) string = "saint_venant_kirchhoff"
end
F = eye(2) + G;
switch material
    case "saint_venant_kirchhoff"
        E = (F.'*F - eye(2))/2;
        P = F*(lambda*trace(E)*eye(2) + 2*mu*E);
    case "neo_hookean"
        FiT = inv(F).';
        P = mu*(F - FiT) + lambda*log(det(F))*FiT;
    otherwise
        error("HyperelasticProblem:material", "Unknown material ""%s"".", material);
end
end
