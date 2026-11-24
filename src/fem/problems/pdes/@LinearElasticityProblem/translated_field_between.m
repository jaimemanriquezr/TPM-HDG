function [V, info] = translated_field_between(obj, old_mesh, new_mesh, degree_from, degree_to, U, shift, options)
% Element-wise L2 projection of a translated piecewise polynomial onto another mesh.
%
%   [V, info] = translated_field_between(obj, old_mesh, new_mesh, degree_from, degree_to, ...
%       U, shift, basis = "dubiner", extension_degree = [])
%
% U holds fields on old_mesh stacked by rows (component blocks), one
% column per old element: Dubiner of degree degree_from (basis =
% "dubiner") or PEERS stresses of degree degree_from = degree_to (basis =
% "peers"). V holds, one column per element K of new_mesh, the dofs of the
% L2 projection onto K's space (Dubiner of degree degree_to, or PEERS of K)
% of x -> U(x - shift). The overlaps K ∩ (K' + shift) with the old
% elements are clipped exactly; on the part of K where x - shift lies in
% no old element, the polynomial of one old element E is extrapolated: the
% old element whose translate overlaps K with the largest area, or, when no
% translate meets K, the one whose translated centroid is nearest to K's
% centroid (lowest index on ties). info.source is E per new element and
% info.uncovered the fraction of K's area where x - shift leaves old_mesh.
% With old_mesh = new_mesh and E = K wherever the fraction is positive
% (small shifts), this is translated_field.
% extension_degree = q in {0, 1} (Dubiner only) replaces, on the
% uncovered part R_K of K, the polynomial of E by its Taylor polynomial of
% degree q about the point of the translate of E nearest to the centroid
% of R_K (info.taylor_point, 2 x N, NaN where K is covered); [] (default)
% or q >= degree_from keeps the polynomial of E.
%
%See also: LinearElasticityProblem.translated_field, LinearElasticityProblem.translate_state
arguments
    obj (1,1) LinearElasticityProblem %#ok<INUSA>
    old_mesh SimplicialGrid
    new_mesh SimplicialGrid
    degree_from (1,1) double
    degree_to (1,1) double
    U double
    shift (2,1) double
    options.basis (1,1) string {mustBeMember(options.basis, ["dubiner", "peers"])} = "dubiner"
    options.extension_degree double {mustBeScalarOrEmpty} = []
end
[V, info] = shifted_projection(old_mesh, degree_from, U, shift, options.basis, new_mesh, degree_to, ...
    options.extension_degree);
end
