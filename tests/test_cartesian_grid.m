%% CartesianGrid: connectivity, geometry and trace parametrization
% On a uniform and a non-uniform grid, checks: entity counts and Euler's
% formula; total area; face lengths; faces run in +x / +y; every face has
% two neighbours (interior) or one duplicated (boundary), and the
% boundary faces are exactly the Dirichlet-tagged ones; the incidence
% table agrees with boundaries; outward normals point away from the
% element centre; the two elements of a face map a shared parameter t to
% the same physical point; and the face-index matrices match the
% numbering. Finally integrates a Q^k function exactly with QkBasis.
clc; clear;
n_fail = 0;
tol = 1e-13;

grids = {CartesianGrid(4, 3), ...
         CartesianGrid([0 0.1 0.35 0.5 1], [-1 -0.2 0.3 2])};
for g = 1:numel(grids)
    mesh = grids{g};
    Nx = mesh.Nx; Ny = mesh.Ny;
    n_el = Nx*Ny; n_f = size(mesh.faces, 2); n_v = size(mesh.coordinates, 2);
    fprintf("grid %d: %d x %d\n", g, Nx, Ny);

    n_fail = n_fail + check("counts", ...
        abs(n_el - size(mesh.elements, 2)) + abs(n_f - (2*Nx*Ny + Nx + Ny)) ...
        + abs(n_v - (Nx+1)*(Ny+1)), 0);
    n_fail = n_fail + check("Euler", abs(n_v - n_f + n_el - 1), 0);

    domain_area = (mesh.x(end) - mesh.x(1)) * (mesh.y(end) - mesh.y(1));
    n_fail = n_fail + check("total area", abs(sum(mesh.volume) - domain_area), tol);

    v1 = mesh.coordinates(:, mesh.faces(1, :));
    v2 = mesh.coordinates(:, mesh.faces(2, :));
    d = v2 - v1;
    n_fail = n_fail + check("face lengths", max(abs(vecnorm(d) - mesh.area)), tol);
    n_fail = n_fail + check("faces +x/+y", ...
        nnz(~((d(1, :) > 0 & d(2, :) == 0) | (d(1, :) == 0 & d(2, :) > 0))), 0);

    count = accumarray(mesh.boundaries(:), 1, [n_f, 1]).';
    is_bdry = count == 1;
    n_fail = n_fail + check("face multiplicity", nnz(count < 1 | count > 2), 0);
    n_fail = n_fail + check("boundary = Dirichlet", ...
        numel(setxor(find(is_bdry), mesh.dirfaces)) + abs(nnz(is_bdry) - 2*(Nx + Ny)), 0);

    % incidence agrees with boundaries
    K = squeeze(mesh.incidence(:, 1, :)); s = squeeze(mesh.incidence(:, 2, :));
    back = mesh.boundaries(sub2ind(size(mesh.boundaries), s, K));
    n_fail = n_fail + check("incidence", nnz(back ~= repmat(1:n_f, 2, 1)), 0);

    % outward normals, and the shared parametrization of each face
    centre = squeeze(mesh.get_physical_points([0.5; 0.5]));
    N = squeeze(mesh.exterior_unit_normals);           % 2 x 4 x n_el
    t = [0.2, 0.7];
    [bad_normal, mismatch] = deal(0);
    for side = 1:4
        rs = mesh.get_2d_local_boundary(t(:), 1, side);
        mid = squeeze(mesh.get_physical_points(mesh.get_2d_local_boundary(0.5, 1, side)));
        bad_normal = bad_normal + nnz(sum(squeeze(N(:, side, :)) .* (mid - centre)) <= 0);
        % physical points of t on the face, from the face's own vertices
        f = mesh.boundaries(side, :);
        on_face = reshape(v1(:, f), 2, 1, []) + reshape(d(:, f), 2, 1, []) .* t;
        mismatch = max(mismatch, max(abs(mesh.get_physical_points(rs) - on_face), [], 'all'));
    end
    n_fail = n_fail + check("outward normals", bad_normal, 0);
    n_fail = n_fail + check("trace param.", mismatch, tol);
    % normals: the face normal [dy; -dx] of each side's face, scaled by
    % its length; outward unit normal = orientation * normal / area
    face_normal = [d(2, :); -d(1, :)];
    [o_err, n_err] = deal(0);
    for side = 1:4
        nf = face_normal(:, mesh.boundaries(side, :));
        n_err = max(n_err, max(abs(mesh.normals(2*side - 1:2*side, :) - nf), [], 'all'));
        o_err = o_err + nnz(any(abs(squeeze(N(:, side, :)) ...
            - mesh.orientation(side, :) .* nf ./ vecnorm(nf)) > tol));
    end
    n_fail = n_fail + check("normals", n_err, tol);
    n_fail = n_fail + check("orientation", o_err, 0);
    n_fail = n_fail + check("orientation S E N W", ...
        nnz(mesh.orientation ~= [1; 1; -1; -1]) + nnz(mesh.permutation ~= 1), 0);

    % face-index matrices
    Fv = mesh.vertical_faces; Fh = mesh.horizontal_faces;
    xv = mesh.coordinates(1, mesh.faces(1, Fv));
    yh = mesh.coordinates(2, mesh.faces(1, Fh));
    n_fail = n_fail + check("face matrices", ...
        nnz(reshape(xv, Nx+1, Ny) ~= mesh.x(:)) + nnz(reshape(yh, Nx, Ny+1) ~= mesh.y), 0);

    % exact integration of u = x^3 y^2 (in Q^3) with QkBasis quadrature
    [x_hat, w] = QkBasis(3).get_quadrature();
    X = mesh.get_physical_points(x_hat);
    u = squeeze(X(1, :, :).^3 .* X(2, :, :).^2);
    integral = sum(w.' * u .* mesh.volume);
    exact = diff(mesh.x([1 end]).^4)/4 * diff(mesh.y([1 end]).^3)/3;
    n_fail = n_fail + check("Q^3 integral", abs(integral - exact), 1e-12);
end
fprintf("%d check(s) failed\n", n_fail);

function failed = check(name, err, tol)
    failed = ~(err <= tol);
    fprintf("  %-22s err = %.2e  %s\n", name, err, string(~failed));
end
