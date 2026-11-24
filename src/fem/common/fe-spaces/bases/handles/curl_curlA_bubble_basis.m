function D_DPB = curl_curlA_bubble_basis(k, r, s, options)
%CURL_CURLA_BUBBLE_BASIS  Bubble modes of the PEERS stress element.
%
%   D_DPB = curl_curlA_bubble_basis(k, r, s, affine_map = G)
%
%   The curl-curl modes built by curl_curl_A_dubiner_bubble_eval from the
%   degree-k Dubiner functions and the element bubble, mapped by the 2-by-2
%   affine_map of each element. For k = 0 there are no modes and a zero
%   column of the right shape is returned. As in dubiner_basis, if s is
%   omitted r is 2-by-n.
    arguments
        k double
        r double
        s (:, 1) double = r(2, :);
        options.affine_map (2, 2, :) double = eye(2);
    end
    if nargin < 3
        r = r(1, :);
    end
    
    r = r(:);
    s = s(:);
    if k == 0
        % no bubble modes: a single zero column keeps the shape, and the
        % element dimension that PEERS_basis reads
        D_DPB = zeros(numel(r), 1, 2, 2, size(options.affine_map, 3));
        return
    end
    D_DPB = curl_curl_A_dubiner_bubble_eval(k, r, s, ...
        inverse_affine_map = options.affine_map);
end
