function spaces = get_discrete_spaces(k)
%-------------------------------------------------------------------------
%Description: FE space descriptors of the Biot problem.
%             Volume:  PEERS P_k            total stress σ;
%                      vector Dubiner P_k^2 displacement η;
%                      scalar Dubiner P_k   rotation γ;
%                      vector Dubiner P_k^2 Darcy flux q;
%                      scalar Dubiner P_k   pore pressure p.
%             Hybrid:  vector Legendre P_k  with 3 components, the trace
%                      (η̂_x, η̂_y, p̂), per face [η̂_x; η̂_y; p̂].
%-------------------------------------------------------------------------
arguments (Input)
    k   double  % polynomial degree
end
arguments (Output)
    spaces  struct  % struct with 'volume' (volume spaces) and 'hybrid' (skeleton spaces)
end
if isempty(k); spaces = struct.empty; return; end
P_k = DubinerBasis(k);
spaces = struct( ...
    "volume", [PEERSBasis(k), P_k.to_vector(), P_k, P_k.to_vector(), P_k], ...
    "hybrid", LegendreBasis(k).to_vector(3));
end
