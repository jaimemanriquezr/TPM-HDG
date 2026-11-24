function mesh = snap_outside(mesh, circle)
%SNAP_OUTSIDE  Interface vertices of a mesh outside a circle moved onto R/cos(a/2).
%
%   mesh = snap_outside(mesh, circle)
%
%   Each vertex of the faces tagged "Interface" (mesh.intfaces) is moved
%   radially to R/cos(a/2), a its widest adjacent interface chord angle,
%   so that every interface chord lies outside the circle, tangent at the
%   midpoint of the widest one: the mesh then lies inside the exterior
%   domain and its transfer paths are O(h^2) and point outward. On a
%   carved mesh this removes the long fluid paths (delta/h_perp up to ~2
%   on a crossed grid) that make the coupled transient unstable at k >= 2.
%
%See also: eulerian_carved_march, ale_translating_march
X = mesh.coordinates;
F = mesh.faces(1:2, mesh.intfaces);
for v = unique(F(:)).'
    a = 0;
    for f = find(any(F == v, 1))
        r = X(:, F(:, f)) - circle.center;
        a = max(a, abs(angle(complex(r(1, 1), r(2, 1))/complex(r(1, 2), r(2, 2)))));
    end
    X(:, v) = circle.center + (X(:, v) - circle.center)/norm(X(:, v) - circle.center)*circle.radius/cos(a/2);
end
mesh.coordinates = X;
x_K = reshape(X(1, mesh.elements), 3, []);
y_K = reshape(X(2, mesh.elements), 3, []);
mesh.volume = .5*sum(x_K.*circshift(y_K, -1) - circshift(x_K, -1).*y_K);
dx_e = diff(reshape(X(1, mesh.faces(1:2, :)), 2, []));
dy_e = diff(reshape(X(2, mesh.faces(1:2, :)), 2, []));
mesh.area = sqrt(dx_e.^2 + dy_e.^2);
n_by_face = [dy_e; -dx_e];
mesh.normals = reshape(n_by_face(:, mesh.boundaries), 2*3, []);
end
