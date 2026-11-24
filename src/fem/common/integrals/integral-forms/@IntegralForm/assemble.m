function [I_A, I_v] = assemble(obj, mesh)
%ASSEMBLE  Integrate an array of integral forms over a mesh.
%
%   [I_A, I_v] = assemble(obj, mesh)
%
%   Sets mesh as the domain of every measure and sums the signed
%   integrals of the terms. I_A collects the terms whose integrand is a
%   BilinearForm and I_v the remaining (linear) terms; a group with no
%   terms is 0. If obj contains no bilinear term both outputs are the
%   assembly of all terms.
is_bilinear_form = arrayfun(@(F) isa(F.integrand, 'BilinearForm'), obj);
if ~any(is_bilinear_form)
    [I_A, I_v] = deal(assemble_form(obj, mesh));
elseif all(is_bilinear_form)
    I_A = assemble_form(obj, mesh);
    I_v = 0;
else
    I_A = assemble_form(obj(is_bilinear_form), mesh);
    I_v = assemble_form(obj(~is_bilinear_form), mesh);
end
end

function I = assemble_form(form, mesh)
if isempty(form)
    I = 0;
    return;
else
    form(1).measure.domain = mesh;
    I = form(1).sign * (form(1).integrand * form(1).measure);
    for i = 2:length(form)
        form(i).measure.domain = mesh;
        I = I + form(i).sign * (form(i).integrand * form(i).measure);
    end
end
end
