function [d, d_prime] = solid_interface_displacement(reference, solid_mesh, dofs, k, options)
%SOLID_INTERFACE_DISPLACEMENT  A solid volume field on Sigma_0, as a function of t.
%
%   [d, d_prime] = solid_interface_displacement(reference, solid_mesh, dofs, k, continuous = false)
%
%   dofs (2*d2 x n_elements, the x block over the y block of Dubiner
%   coefficients, as state.solid.displacement) is evaluated at X(t) on the
%   reference curve, each t on the solid interface face whose segment
%   (map_mesh_to_interface) contains it, through the polynomial of that
%   face's element extrapolated to X(t): the TPM evaluation of the solid
%   field on Sigma_0. d(t) and d_prime(t) take a row of parameters and
%   return 2 x n; d_prime is the fourth-order central difference of the
%   same element's polynomial (the piece of t is kept across the stencil),
%   so d is piecewise smooth and d_prime exact up to O(h^4) on each piece.
%   For DeformedCurve(reference, d, parametric = true, derivative = d_prime).
%
%   The pieces of neighbouring elements disagree at their common
%   breakpoint by the jump of the discrete field, so X + d is a
%   discontinuous curve. continuous = true adds to each piece the linear
%   function of t that takes it, at both ends, to the average of the two
%   one-sided values there (the ends of an open curve are left as they
%   are): d is then continuous, changed by at most the largest jump, and
%   d_prime gains the slope of that correction.
%
%See also: DeformedCurve, lagrangian_carved_march
arguments
    reference
    solid_mesh
    dofs
    k (1,1) double
    options.continuous (1,1) logical = false
end
segments = reference.map_mesh_to_interface(solid_mesh);
[t_lo, order] = sort(segments.t(1, :));
t_hi = segments.t(2, order);
faces = segments.face(order);
K = solid_mesh.get_element_by_face(faces);
K = K(1, :);
t0 = reference.t_range(1);
period = diff(reference.t_range);
periodic = is_periodic(reference);
if periodic
    t_hi(t_hi < t_lo) = t_hi(t_hi < t_lo) + period;
end
d2 = nchoosek(k + 2, 2);
h = 2e-4*period;
n_pieces = numel(t_lo);

% the linear correction of each piece: c_lo at t_lo, c_hi at t_hi
[c_lo, c_hi] = deal(zeros(2, n_pieces));
if options.continuous && n_pieces > 1
    left = evaluate(t_hi, 1:n_pieces);       % piece e at its right end
    right = evaluate(t_lo, 1:n_pieces);      % piece e at its left end
    following = [2:n_pieces, 1];
    mean_value = (left + right(:, following))/2;
    c_hi = mean_value - left;
    c_lo(:, following) = mean_value - right(:, following);
    if ~periodic
        c_hi(:, end) = 0;
        c_lo(:, 1) = 0;
    end
end
slope = (c_hi - c_lo)./(t_hi - t_lo);

d = @(t) evaluate(t, piece(t)) + correction(t, piece(t));
d_prime = @(t) (evaluate(t - 2*h, piece(t)) - 8*evaluate(t - h, piece(t)) ...
    + 8*evaluate(t + h, piece(t)) - evaluate(t + 2*h, piece(t))) / (12*h) + slope(:, piece(t));

    function e = piece(t)
        if periodic, t = t0 + mod(t - t0, period); end
        e = discretize(reshape(t, 1, []), [-Inf, t_lo(2:end), Inf]);
    end

    function c = correction(t, e)
        t = reshape(t, 1, []);
        if periodic, t = t0 + mod(t - t0, period); end
        c = c_lo(:, e) + slope(:, e).*(t - t_lo(e));
    end

    function v = evaluate(t, e)
        t = reshape(t, 1, []);
        X = reference.get_2d_points(t);
        rs = solid_mesh.get_2d_local(reshape(X, 2, 1, []), K(e).');
        B = dubiner_basis(k, reshape(rs(1, 1, :), [], 1), reshape(rs(2, 1, :), [], 1));
        v = [sum(B .* dofs(1:d2, K(e)).', 2), sum(B .* dofs(d2 + (1:d2), K(e)).', 2)].';
    end
end
