%-------------------------------------------------------------------------
%Description: Signed TPM sliver integral against the mesh's area deficit.
%
%             Contracted with the patch/path quadrature weights, the
%             get_sliver_jacobian area element must give the area the
%             mesh's computational domain is missing (positive) or covers
%             in excess (negative): deficit = |subdomain| - sum(volume),
%             for either mesh of the pair and either side of Sigma.
%             "jacobian" is the exact signed area up to the patch
%             quadrature (observed relative error <= 3.3e-12, checked at
%             1e-10); "rectangular" approximates it, so only its sign is
%             checked. The kidney (a LevelSetCurve) checks the jacobian
%             against a shoelace area of Sigma; its parametrization is
%             arc length only up to the contour discretization (observed
%             relative error <= 5e-8), hence the looser 1e-6.
%-------------------------------------------------------------------------
n_failed = 0;

sq = 4;  tq = 5;
[~, wS] = lgwt(sq, 0, 1);
[~, wTP] = lgwt(tq, 0, 1);
integrate = @(curve, mesh, part, pw) sum((wS(:) .* get_sliver_jacobian(curve, ...
    mesh, part, patch_weight = pw, quadrature_patch = sq, quadrature_path = tq)) ...
    .* wTP(:)', 'all');

% circle: exact subdomain areas
r = .25;
curve = CircularCurve([.5; .5], r);
for family = ["interior", "fitted", "protrude_fluid"]
    for side = ["outer", "inner"]
        [fm, sm] = interface_meshes(curve, family, 0.05, fluid_side = side);
        [pF, pS] = get_partitions(curve, fm, sm, quadrature_patch = sq, quadrature_path = tq);
        A_inner = pi * r^2;
        A_fluid = A_inner + (side == "outer") * (1 - 2 * A_inner);
        meshes = {"fluid", fm, pF, A_fluid; "solid", sm, pS, 1 - A_fluid};
        for i = 1:2
            [name, mesh, part, A] = meshes{i, :};
            deficit = A - sum(mesh.volume);
            label = sprintf("circle %s/%s, %s", family, side, name);
            I = integrate(curve, mesh, part, "jacobian");
            fprintf("  %-40s deficit %+.6e  jacobian %+.6e  rectangular %+.6e\n", ...
                label, deficit, I, integrate(curve, mesh, part, "rectangular"));
            n_failed = n_failed + check(label + ": jacobian = deficit", ...
                abs(I - deficit) <= 1e-10 * abs(deficit));
            n_failed = n_failed + check(label + ": sign(rectangular) = sign(deficit)", ...
                sign(integrate(curve, mesh, part, "rectangular")) == sign(deficit));
        end
    end
end

% kidney (level-set curve): area of Sigma from a fine shoelace
try
    curve = discretize(KidneyCurve(), n_grid = 4096);
    P = curve.get_2d_points(linspace(curve.t_range(1), curve.t_range(2), 2^16 + 1));
    P = P(:, 1:end-1);
    A_kidney = abs(sum(P(1, :) .* circshift(P(2, :), -1) - circshift(P(1, :), -1) .* P(2, :))) / 2;
    [fm, sm] = interface_meshes(curve, "interior", 1/12, box = [-1 1 -1 1], fluid_side = "outer");
    [pF, pS] = get_partitions(curve, fm, sm, quadrature_patch = sq, quadrature_path = tq);
    meshes = {"fluid", fm, pF, 4 - A_kidney; "solid", sm, pS, A_kidney};
    for i = 1:2
        [name, mesh, part, A] = meshes{i, :};
        deficit = A - sum(mesh.volume);
        label = "kidney interior/outer, " + name;
        I = integrate(curve, mesh, part, "jacobian");
        fprintf("  %-40s deficit %+.6e  jacobian %+.6e\n", label, deficit, I);
        n_failed = n_failed + check(label + ": jacobian = deficit", ...
            abs(I - deficit) <= 1e-6 * abs(deficit));
    end
catch err
    fprintf(2, "  kidney: %s\n", err.message);
    n_failed = n_failed + check("kidney: runs", false);
end

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
