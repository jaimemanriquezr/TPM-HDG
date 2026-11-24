function varargout = get_test_functions(obj)
% Test functions of all the problem's spaces, each carrying the
% problem's mesh (see get_trial_functions).
[varargout{1:nargout}] = get_test_functions(obj.Space);
for i = 1:numel(varargout)
    varargout{i}.mesh = obj.Mesh;
end
end
