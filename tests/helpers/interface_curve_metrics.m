function m = interface_curve_metrics(position, breaks, exact_position, geometry, options)
%INTERFACE_CURVE_METRICS  Geometric errors of a piecewise reconstructed interface.
%
%   m = interface_curve_metrics(position, breaks, exact_position, geometry, ...
%       n_per_piece = 9, period = 1)
%
%   position        handle s -> 2 x n, the reconstructed curve on the
%                   material parameter, e.g. @(s) reference.get_2d_points(s) + D(s).
%   breaks          the piece breaks (info.breaks of
%                   eulerian_interface_displacement): 1 x (n_pieces + 1) for
%                   an open curve, the piece starts for a closed one (the
%                   last piece ends at breaks(1) + period).
%   exact_position  handle s -> 2 x n, the exact material point x(s), or [].
%   geometry        the exact curve as a set: struct("type", "flat",
%                   "height", y) (open) or struct("type", "circle",
%                   "center", c, "radius", r) (closed).
%
%   Each piece [a, b] is sampled at the n_per_piece Chebyshev-Lobatto
%   points of [a + delta, b - delta], delta = 1e-12 (b - a), so that the
%   samples of a piece, its ends included, come from its own polynomial.
%
%   m.distance        max over the samples of the distance to the exact
%                     set (flat |y - height|, circle ||x - c| - r|).
%   m.point_error     max |position - exact_position| (NaN without one).
%   m.measure         flat: the mean height, int y dx / (x_end - x_start);
%                     circle: the enclosed area, (1/2) int (x dy - y dx)
%                     about the centre.
%   m.measure_ratio   mean height / height, or area / (pi r^2).
%   m.monotonicity    the smallest increment, over consecutive samples
%                     (across breaks too, and across the wrap), of the
%                     advance along the exact curve (flat: x; circle: the
%                     unwrapped polar angle about the centre), over the
%                     mean increment; > 0 when the curve is traversed
%                     monotonically. An increment across a break is about
%                     2 delta |x'| plus the jump's advance, so the value is
%                     tiny but positive for a curve without fold-back.
%   m.max_jump        the largest |jump| of position at the interior
%                     breaks (closed: the wrap included).
%
%   Quadrature. On each piece the curve is smooth, and the integrals of
%   y x' (flat) and x y' - y x' (circle) are computed by Clenshaw-Curtis
%   on the samples, with x' from the Chebyshev differentiation matrix of
%   the same samples: no tangent is needed, and the rule is spectrally
%   accurate in n_per_piece on each piece (exact for the polynomial
%   interpolant of degree n_per_piece - 1; on a piece of length h the
%   error is O(h^n_per_piece)). The gaps between pieces (the jumps and the
%   2 delta of the sampling) are closed by straight chords (trapezoid /
%   triangle), so a discontinuous curve encloses the polygon-closed area.
%   On an exact smooth curve with 9 points per piece and pieces of a
%   tenth of the parameter range the relative error is below 1e-10.
%
%See also: eulerian_interface_displacement
arguments
    position function_handle
    breaks (1,:) double
    exact_position
    geometry struct
    options.n_per_piece (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.n_per_piece, 3)} = 9
    options.period (1,1) double = 1
end
closed = geometry.type == "circle";
if closed
    s_lo = breaks;
    s_hi = [breaks(2:end), breaks(1) + options.period];
else
    s_lo = breaks(1:end - 1);
    s_hi = breaks(2:end);
end
n_pieces = numel(s_lo);
n = options.n_per_piece;
[nodes, Dm, w] = chebyshev(n - 1);
delta = 1e-12*(s_hi - s_lo);
a = s_lo + delta;
b = s_hi - delta;
s = a + (nodes + 1)/2 .* (b - a);              % n x n_pieces
P = reshape(position(reshape(s, 1, [])), 2, n, n_pieces);

% distance to the exact set and to the exact points
if closed
    xi = P - geometry.center(:);
    distance = abs(vecnorm(xi, 2, 1) - geometry.radius);
else
    distance = abs(P(2, :, :) - geometry.height);
end
m.distance = max(distance, [], "all");
m.point_error = NaN;
if ~isempty(exact_position)
    E = reshape(exact_position(reshape(s, 1, [])), 2, n, n_pieces);
    m.point_error = max(vecnorm(P - E, 2, 1), [], "all");
end

% one-sided ends and the chords closing the gaps between pieces
first = squeeze(P(:, 1, :));
last = squeeze(P(:, n, :));
if closed
    next = [2:n_pieces, 1];
else
    next = 2:n_pieces;
end
from = last(:, 1:numel(next));
to = first(:, next);
m.max_jump = max([0, vecnorm(to - from, 2, 1)]);

% the measure: Clenshaw-Curtis per piece, chords across the gaps
x1 = squeeze(P(1, :, :));
x2 = squeeze(P(2, :, :));
if closed
    c = geometry.center(:);
    x1 = x1 - c(1);  x2 = x2 - c(2);
    pieces = w*(x1.*(Dm*x2) - x2.*(Dm*x1));
    chords = (from(1, :) - c(1)).*(to(2, :) - c(2)) - (from(2, :) - c(2)).*(to(1, :) - c(1));
    m.measure = (sum(pieces) + sum(chords))/2;
    m.measure_ratio = m.measure/(pi*geometry.radius^2);
else
    pieces = w*(x2.*(Dm*x1));
    chords = (from(2, :) + to(2, :))/2 .* (to(1, :) - from(1, :));
    m.measure = (sum(pieces) + sum(chords))/(last(1, end) - first(1, 1));
    m.measure_ratio = m.measure/geometry.height;
end

% monotone advance along the exact curve
if closed
    advance = unwrap(atan2(x2(:), x1(:)));
    advance(end + 1) = advance(1) + 2*pi;
else
    advance = x1(:);
end
increment = diff(advance);
m.monotonicity = min(increment)/mean(increment);
end

function [x, Dm, w] = chebyshev(N)
% Chebyshev-Lobatto points of [-1, 1] in increasing order, their
% differentiation matrix and the Clenshaw-Curtis weights (Trefethen,
% Spectral Methods in MATLAB, cheb and clencurt)
theta = pi*(N:-1:0).'/N;
x = cos(theta);
c = [2; ones(N - 1, 1); 2].*(-1).^(0:N).';
dX = x - x.';
Dm = (c*(1./c).')./(dX + eye(N + 1));
Dm = Dm - diag(sum(Dm, 2));
w = zeros(1, N + 1);
inner = 2:N;
v = ones(N - 1, 1);
if mod(N, 2) == 0
    w([1, N + 1]) = 1/(N^2 - 1);
    for j = 1:N/2 - 1
        v = v - 2*cos(2*j*theta(inner))/(4*j^2 - 1);
    end
    v = v - cos(N*theta(inner))/(N^2 - 1);
else
    w([1, N + 1]) = 1/N^2;
    for j = 1:(N - 1)/2
        v = v - 2*cos(2*j*theta(inner))/(4*j^2 - 1);
    end
end
w(inner) = 2*v/N;
end
