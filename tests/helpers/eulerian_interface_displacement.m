function [D, D_prime, info] = eulerian_interface_displacement(reference, curve, solid_mesh, dofs, k, rigid, options)
%EULERIAN_INTERFACE_DISPLACEMENT  The interface displacement of an Eulerian solid, on Sigma_0.
%
%   [D, D_prime, info] = eulerian_interface_displacement(reference, curve, solid_mesh, ...
%       dofs, k, rigid, continuous = false)
%
%   The solid occupies the current carve solid_mesh; its total Eulerian
%   displacement is w(x) = rigid + u(x), u the elastic displacement with
%   Dubiner dofs of degree k (2*d2 x n_elements, the x block over the y
%   block), and the particle at x comes from x - w(x). For the material
%   parameter s of the reference curve Sigma_0, X(s), the current position
%   x_h(s) solves
%       x - rigid - u_K(x) = X(s),
%   u_K the polynomial of one element K extrapolated to x. D(s) = x_h(s) -
%   X(s) and D_prime(s) = x_h'(s) - X'(s), x_h' = (I - grad u_K(x_h))^-1
%   X'(s) (the Dubiner gradients through the affine map), take a row of
%   parameters and return 2 x n; for
%   DeformedCurve(reference, D, parametric = true, derivative = D_prime).
%
%   Pieces. curve is the curve the interface faces of solid_mesh are
%   mapped to (map_mesh_to_interface), with a parameter sigma that need
%   not be material but agrees with s to within half a period on a closed
%   curve. Each solid interface face e covers [sigma_lo, sigma_hi] and has
%   its element K(e) (a face split at the seam of a closed curve is one
%   piece, continued past the seam). At the point p_j of curve between two consecutive
%   faces, each of the two elements gives a material parameter, the
%   closest point on the reference of p_j - rigid - u_K(p_j); their mean
%   is the break s_j, and the piece [s_j, s_j+1] is solved with K(e)
%   throughout (slightly extrapolated along the curve near its ends). The
%   pieces are chosen where x_h is, not where X is, since u may slide the
%   material tangentially by more than an element. The ends of an open
%   curve are those of reference.t_range. Breaks that do not increase
%   raise eulerian_interface_displacement:foldedPieces.
%
%   x_h is found by Newton on F(x) = x - rigid - u_K(x) - X(s), from the
%   point of curve at the sigma linearly interpolated from the piece's
%   ends, to an absolute step below 1e-14 (at most 30 iterations, else
%   eulerian_interface_displacement:notConverged).
%
%   Neighbouring pieces disagree at their break by the jump of the
%   discrete field. continuous = true adds to each piece the linear
%   function of s that takes it, at both ends, to the average of the two
%   one-sided values there (the ends of an open curve are left as they
%   are; a closed curve wraps), as solid_interface_displacement does; its
%   slope is added to D_prime.
%
%   info.breaks    the s_j: 1 x (n_pieces + 1) on an open curve, the piece
%                  starts (1 x n_pieces, the last piece ending at
%                  breaks(1) + period) on a closed one.
%   info.jumps     2 x n_breaks, the jump (right minus left) of the
%                  uncorrected x_h at each interior break (closed: at
%                  every piece start, the first across the wrap).
%   info.elements  K(e), 1 x n_pieces.
%   info.on_piece  handle (s, e) -> D(s) evaluated with piece e (e scalar
%                  or one per parameter), correction included: the
%                  one-sided values at a break.
%
%   The Dubiner evaluation is not reliable at an extrapolated point lying
%   on the reference line s = 1 through the collapsed vertex.
%
%See also: solid_interface_displacement, DeformedCurve, eulerian_solid_march
arguments
    reference
    curve
    solid_mesh
    dofs
    k (1,1) double
    rigid (2,1) double
    options.continuous (1,1) logical = false
end
segments = curve.map_mesh_to_interface(solid_mesh);
[sigma_lo, order] = sort(segments.t(1, :));
sigma_hi = segments.t(2, order);
faces = segments.face(order);
K = solid_mesh.get_element_by_face(faces);
K = K(1, :);
periodic = is_periodic(reference);
period = diff(reference.t_range);
if periodic
    sigma_hi(sigma_hi < sigma_lo) = sigma_hi(sigma_hi < sigma_lo) + period;
    if numel(faces) > 1 && faces(1) == faces(end)
        % a face split at the seam is one piece, continued past the seam
        sigma_hi(end) = sigma_hi(1) + period;
        [sigma_lo, sigma_hi, K] = deal(sigma_lo(2:end), sigma_hi(2:end), K(2:end));
    end
end
n_pieces = numel(sigma_lo);
d2 = nchoosek(k + 2, 2);

% per piece: the inverse affine map, the image of the reference origin
% and the coefficients (n_pieces x d2 per component)
G = get_affine_map(solid_mesh, K, inverse = true);
G = reshape(G, 4, n_pieces);                     % G11 G21 G12 G22
[~, x3] = get_affine_map(solid_mesh, K);
C1 = dofs(1:d2, K).';
C2 = dofs(d2 + (1:d2), K).';

% the breaks: the material parameter of each piece start, from either side
if periodic
    starts = 1:n_pieces;
    left_of = [n_pieces, 1:n_pieces - 1];
else
    starts = 2:n_pieces;
    left_of = 1:n_pieces - 1;
end
p = curve.get_2d_points(sigma_lo(starts));
s_left = material(p - rigid - field(p, left_of));
s_right = material(p - rigid - field(p, starts));
if periodic
    % the two one-sided parameters on one branch (a break at the seam may
    % have them on either side of it)
    wrapped = abs(s_right - s_left) > period/2;
    s_right(wrapped) = s_right(wrapped) - period*round((s_right(wrapped) - s_left(wrapped))/period);
end
s_break = (s_left + s_right)/2;
if periodic
    % the branch nearest the curve parameter
    s_break = s_break + period*round((sigma_lo(starts) - s_break)/period);
    breaks = s_break;
    s_end = [breaks(2:end), breaks(1) + period];
else
    breaks = [reference.t_range(1), s_break, reference.t_range(2)];
    s_end = breaks(2:end);
end
s_start = breaks(1:n_pieces);
if any(s_end - s_start <= 0)
    error("eulerian_interface_displacement:foldedPieces", ...
        "The material breaks of the pieces do not increase (%d of %d pieces folded).", ...
        nnz(s_end - s_start <= 0), n_pieces);
end

% one-sided values, jumps and the linear correction of each piece
left = position(s_end, 1:n_pieces);          % piece e at its right end
right = position(s_start, 1:n_pieces);       % piece e at its left end
following = [2:n_pieces, 1];
if periodic
    jumps = right - left(:, left_of);
else
    jumps = right(:, 2:end) - left(:, 1:end - 1);
end
[c_lo, c_hi] = deal(zeros(2, n_pieces));
if options.continuous && n_pieces > 1
    mean_value = (left + right(:, following))/2;
    c_hi = mean_value - left;
    c_lo(:, following) = mean_value - right(:, following);
    if ~periodic
        c_hi(:, end) = 0;
        c_lo(:, 1) = 0;
    end
end
slope = (c_hi - c_lo)./(s_end - s_start);

D = @(s) displacement(s, []);
D_prime = @(s) derivative(s);
info = struct("breaks", breaks, "jumps", jumps, "elements", K, ...
    "on_piece", @(s, e) displacement(s, e));

    function s = wrap(s)
        s = reshape(s, 1, []);
        if periodic, s = breaks(1) + mod(s - breaks(1), period); end
    end

    function e = piece(s)
        e = discretize(s, [-Inf, breaks(2:n_pieces), Inf]);
    end

    function v = displacement(s, e)
        s = wrap(s);
        if isempty(e), e = piece(s); end
        e = e + 0*s;
        if periodic
            % the branch of s nearest piece e (its right end, not the wrap)
            s = s + period*round(((s_start(e) + s_end(e))/2 - s)/period);
        end
        v = position(s, e) - reference.get_2d_points(s) + c_lo(:, e) + slope(:, e).*(s - s_start(e));
    end

    function v = derivative(s)
        s = wrap(s);
        e = piece(s);
        [~, J] = position(s, e);
        T = reference.tangent(s) + 0*s;
        v = solve2(J, T) - T + slope(:, e);
    end

    function [x, J] = position(s, e)
        % x_h(s) on piece e by Newton, and I - grad u_K there (4 x n)
        X = reference.get_2d_points(s);
        sigma = sigma_lo(e) + (s - s_start(e))./(s_end(e) - s_start(e)).*(sigma_hi(e) - sigma_lo(e));
        x = reshape(curve.get_2d_points(sigma), 2, []);
        for it = 1:30
            [u, J] = field(x, e);
            step = solve2(J, X + rigid + u - x);
            x = x + step;
            if max(abs(step), [], "all") < 1e-14
                if nargout > 1, [~, J] = field(x, e); end
                return
            end
        end
        error("eulerian_interface_displacement:notConverged", ...
            "Newton did not reach a step below 1e-14 in 30 iterations (last %.2e).", ...
            max(abs(step), [], "all"));
    end

    function [u, J] = field(x, e)
        % u_K(e) at the columns of x, and I - grad u_K (rows J11 J21 J12 J22)
        y = x - x3(:, e);
        r = (G(1, e).*y(1, :) + G(3, e).*y(2, :)).';
        q = (G(2, e).*y(1, :) + G(4, e).*y(2, :)).';
        if nargout > 1
            [B, Br, Bq] = dubiner_basis(k, r, q);
        else
            B = dubiner_basis(k, r, q);
        end
        a = C1(e, :);
        b = C2(e, :);
        u = [sum(B.*a, 2), sum(B.*b, 2)].';
        if nargout > 1
            ar = sum(Br.*a, 2).';  aq = sum(Bq.*a, 2).';
            br = sum(Br.*b, 2).';  bq = sum(Bq.*b, 2).';
            J = [1 - (ar.*G(1, e) + aq.*G(2, e)); -(br.*G(1, e) + bq.*G(2, e)); ...
                -(ar.*G(3, e) + aq.*G(4, e)); 1 - (br.*G(3, e) + bq.*G(4, e))];
        end
    end

    function s = material(q)
        % parameter of the closest point of the reference, polished by
        % Gauss-Newton steps
        s = zeros(1, size(q, 2));
        for i = 1:numel(s)
            s(i) = reference.find_closest_point(q(1, i), q(2, i));
        end
        for it = 1:3
            T = reference.tangent(s) + 0*s;
            s = s - sum((reference.get_2d_points(s) - q).*T, 1)./sum(T.^2, 1);
        end
    end
end

function z = solve2(J, b)
% J \ b columnwise, J stored as rows J11 J21 J12 J22
det = J(1, :).*J(4, :) - J(3, :).*J(2, :);
z = [J(4, :).*b(1, :) - J(3, :).*b(2, :); J(1, :).*b(2, :) - J(2, :).*b(1, :)]./det;
end
