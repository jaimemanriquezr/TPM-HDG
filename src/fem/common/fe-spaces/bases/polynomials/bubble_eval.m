function [p, dpr, dps] = bubble_eval(r, s)
    % Normalized cubic bubble b = 27 lambda_1 lambda_2 lambda_3
    % = 27 r s (1 - r - s), value 1 at the barycenter. The toolbox-wide
    % convention: curl_curlA_bubble_basis carries the same factor, so both
    % PEERS bases produce identical coefficients.
    p = 27 * r .* s .* (1 - r - s);
    dpr = 27 * (s .* (1 - 2*r - s));
    dps = 27 * (r .* (1 - r - 2*s));
end
