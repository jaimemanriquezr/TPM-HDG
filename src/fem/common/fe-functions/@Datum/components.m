function comps = components(obj, n_comp)
    % split a vector-valued Datum into n_comp scalar component Datums
    % (a cell column); a constant is repeated in every component
    v = obj.value;
    switch class(v)
        case 'cell'
            assert(numel(v) == n_comp, ...
                'Datum has %d components, trace space expects %d.', numel(v), n_comp);
            comps = cellfun(@(f) {Datum(f)}, v(:));
        case 'function_handle'
            assert(n_comp == 1);
            comps = {obj};
        case 'double'
            comps = arrayfun(@(a) {Datum(a)}, repmat(v(:), n_comp, 1));
            comps = comps(1:n_comp);
        otherwise
            error('Datum:components', ...
                  'Unsupported datum value of class %s.', class(v));
    end
end
