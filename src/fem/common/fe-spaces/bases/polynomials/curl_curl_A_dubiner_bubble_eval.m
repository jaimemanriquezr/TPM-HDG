function D_DPB = curl_curl_A_dubiner_bubble_eval(k, r, s, options)
%-------------------------------------------------------------------------
%Description: Evaluates the k + 1 matrix-valued bubble modes of the PEERS
%           basis at the natural coordinates (r, s), built from the first
%           and second derivatives of the k + 1 highest-index Dubiner
%           polynomials and the cubic bubble of bubble_eval. Returns []
%           for k = 0.
%Input:     - k: polynomial order
%           - r, s: column vectors of natural coordinates
%           - inverse_affine_map: 2-by-2-by-n_pages inverse affine maps
%             that carry the derivatives to the physical element
%             (default eye(2))
%Output:    - D_DPB: values, N-by-(k+1)-by-2-by-2-by-n_pages
%-------------------------------------------------------------------------
    arguments 
        k double
        r (:, 1) double
        s (:, 1) double = r(2, :);
        options.inverse_affine_map (2, 2, :) double = eye(2);
    end
    if nargin < 3; r = r(1, :); end
    if k == 0; D_DPB = []; return; end
    
    r = r(:);
    s = s(:);

    d2 = nchoosek(k+2, 2);
    [~, d_Px, d_Py, dd_Pxx, dd_Pyy, dd_Pxy] = DubinerBasis.evaluate_basis((d2-k):d2, r, s);
    [B, DBx, DBy] = bubble_eval(r, s);
    DB_K = [DBx, DBy];
    
    G = options.inverse_affine_map;
    DP = deal(zeros([size(d_Px) 2 size(G, 3)]));
    DB = deal(zeros([size(DB_K) size(G, 3)]));
    [D_DPB, DDP] = deal(zeros([size(dd_Pxx) 2 2 size(G, 3)]));
    for i = 1:2
        DP(:, :, i, :) = G(1,i,:).*d_Px + G(2,i,:).*d_Py;
        DB(:, i, :) = G(1,i,:).*DB_K(:,1) + G(2,i,:).*DB_K(:,2);
        for j = 1:2
            DDP(:, :, i, j, :) = dd_Pxx .* (G(1,i,:).*G(1,j,:)) ...
                                + dd_Pxy .* (G(2,i,:).*G(1,j,:) + G(1,i,:).*G(2,j,:)) ...
                                + dd_Pyy .* (G(2,i,:).*G(2,j,:));
        end
    end
    DB = reshape(DB, [size(DB, 1) 1 size(DB, [2 3])]);

    s_plus = @(A, B) squeeze(A) + squeeze(B);
    D_DPB(:, :, 1, 1, :) = s_plus(DDP(:,:,1,2,:).*B, DP(:,:,1,:).*DB(:,:,2,:));
    D_DPB(:, :, 2, 1, :) = s_plus(DDP(:,:,2,2,:).*B, DP(:,:,2,:).*DB(:,:,2,:));
    D_DPB(:, :, 1, 2, :) = -s_plus(DDP(:,:,1,1,:).*B, DP(:,:,1,:).*DB(:,:,1,:));
    D_DPB(:, :, 2, 2, :) = -s_plus(DDP(:,:,2,1,:).*B, DP(:,:,2,:).*DB(:,:,1,:));
end