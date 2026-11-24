%-------------------------------------------------------------------------
%Description: FEField.compute_l2_error on tensor, vector and scalar fields.
%
%             The L2 norm of a tensor field is the Frobenius norm of its
%             components, so a zero field against the constant a*I, or
%             against the skew [0 a; -a 0], on the unit square must give
%             a*sqrt(2). The vector and scalar cases pin the Euclidean
%             combination of components.
%-------------------------------------------------------------------------
n_failed = 0;

mesh = StructuredSimplicialGrid(2, 2);
n_el = size(mesh.elements, 2);
a = 0.7;
zero = @(x, y) 0*x;
const = @(c) @(x, y) c + 0*x;

sigma = FEField(PEERSBasis(1), zeros(14, n_el));
u = FEField(DubinerBasis(1).to_vector(), zeros(6, n_el));
p = FEField(DubinerBasis(1), zeros(3, n_el));

e = compute_l2_error(sigma, mesh, Datum({const(a), zero; zero, const(a)}));
n_failed = n_failed + check("tensor a*I", abs(e - a*sqrt(2)) < 1e-14);
e = compute_l2_error(sigma, mesh, Datum({zero, const(a); const(-a), zero}));
n_failed = n_failed + check("tensor skew", abs(e - a*sqrt(2)) < 1e-14);
e = compute_l2_error(u, mesh, Datum({const(3); const(4)}));
n_failed = n_failed + check("vector (3, 4)", abs(e - 5) < 1e-14);
e = compute_l2_error(p, mesh, Datum(const(a)));
n_failed = n_failed + check("scalar", abs(e - a) < 1e-14);

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
