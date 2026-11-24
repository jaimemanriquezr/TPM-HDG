function V = translated_field(~, mesh, degree, U, shift)
% Element-wise L2 projection of a translated piecewise polynomial.
%
%   V = translated_field(obj, mesh, degree, U, shift)
%
% U holds Dubiner fields of the given degree stacked by rows (component
% blocks), one column per element of mesh; V holds the dofs of the L2
% projection onto the same space of x -> U(x - shift), the translate
% translate_state applies to the volume levels (shifted_projection: exact
% overlaps, K's own polynomial extrapolated where x - shift leaves the
% mesh).
%
%See also: LinearElasticityProblem.translate_state
V = shifted_projection(mesh, degree, U, shift, "dubiner");
end
