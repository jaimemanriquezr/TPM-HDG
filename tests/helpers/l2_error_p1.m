function e = l2_error_p1(mesh, k, u_star, u_exact)
% L2 error of a P_{k+1} Dubiner vector field against handles
[xy, rs, w] = quadpoints_elements_2d(mesh, k + 2, output = "cell");
P1 = dubiner_basis(k + 1, rs{:});
e = 0;
for i = 1:2
    exact = u_exact{i}(xy{:}) + 0*xy{1};
    e = e + sum(mesh.volume .* sum(w .* (P1*u_star{i} - exact).^2, 1));
end
e = sqrt(e);
end
