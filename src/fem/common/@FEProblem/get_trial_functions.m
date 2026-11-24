function varargout = get_trial_functions(obj)
% Trial functions of all the problem's spaces, each carrying the
% problem's mesh -- required by the transmission forms, which route
% every term by its trial's own mesh (Measure(mesh, "dI") design).
[varargout{1:nargout}] = get_trial_functions(obj.Space);
for i = 1:numel(varargout)
    varargout{i}.mesh = obj.Mesh;
end
end
