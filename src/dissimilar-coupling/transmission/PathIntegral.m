function obj = PathIntegral(x, options)
% PathIntegral(u): the transfer-path integral
%   int_0^1 u(x + t m) . m dt
% of a local trial expression along its side's transfer paths; constant
% matrix weights spell directly on the argument (PathIntegral(kappa\u)).
% Returns a ExtendedFunction (see Extrapolate on why a function).
arguments
    x
    options.weight double = []   % deprecated; prefer weights on x
end
obj = ExtendedFunction(x);
for i = 1:numel(obj.terms)
    assert(obj.terms(i).kind == "onface", "PathIntegral:arg", ...
        "PathIntegral wraps local (volume) trials only");
    obj.terms(i).kind = "path";
    if ~isempty(options.weight)
        assert(isempty(obj.terms(i).weight), "PathIntegral:weight", ...
            "term already carries a weight; do not also pass weight=");
        obj.terms(i).weight = options.weight;
    end
end
end
