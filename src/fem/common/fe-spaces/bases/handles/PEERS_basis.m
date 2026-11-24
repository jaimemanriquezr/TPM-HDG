function P = PEERS_basis(k, r, s, G)
%PEERS_BASIS  PEERS stress basis of degree k on the reference triangle.
%
%   P = PEERS_basis(k, r, s, G)
%
%   The scalar Dubiner basis times each of the four elementary 2-by-2
%   matrices, followed (for k > 0) by the bubble modes of
%   curl_curlA_bubble_basis, with G the 2-by-2 affine map of each element.
%   If s is omitted r is 2-by-n.
    arguments
        k double
        r double
        s double = r(2, :);
        G (2, 2, :) double = eye(2);
    end
    if nargin < 3
        r = r(1, :);
    end

    p = dubiner_basis(k, r, s);
    pp = [p.*reshape([1, 0; 0, 0], [1 1 2 2]), ...
          p.*reshape([0, 0; 1, 0], [1 1 2 2]), ...
          p.*reshape([0, 1; 0, 0], [1 1 2 2]), ...
          p.*reshape([0, 0; 0, 1], [1 1 2 2])];

    ccAb = curl_curlA_bubble_basis(k, r, s, affine_map = G);
    pp_rep = repmat(pp, [1 1 1 1 size(ccAb, 5)]);
    if k > 0
        P = cat(2, pp_rep, ccAb);
    else
        P = pp_rep;  % bubble degenerates for k=0 (grad P_0 = 0)
    end
end