function [w, mu, varargout] = evaluate_quadrature(mesh, u, options)
    %EVALUATE_QUADRATURE  Quadrature weights, measures and factor values on a mesh.
    %
    %   [w, mu, u1, u2, ...] = evaluate_quadrature(mesh, f1, f2, ..., ...
    %       measure_type = 'dx', weight = 1, order = [], faces = [])
    %
    %   evaluates the factors f1, f2, ... (FEFields or composite nodes) at
    %   the reference quadrature points of the integration domain and
    %   returns the reference weights w, the measure mu of each entity
    %   multiplied by weight, and the evaluations u1, u2, ... of the
    %   factors. measure_type selects the domain: 'dx' elements (mu =
    %   element areas), 'dS' the sides of every element (mu = the
    %   face lengths of mesh.boundaries), 'ds' mesh faces, restricted to
    %   faces when given, with both adjacent elements evaluated (last
    %   dimension 2), and 'de' the factors evaluated with the 1-D rule
    %   (mu = face lengths). order defaults to the largest order among
    %   the factors.
    %
    %See also: AbstractFEFunction.integrate, get_1d_reference_quadrature, get_2d_reference_quadrature
    arguments
        mesh AbstractGrid = SimplicialGrid(reference_element(), ...
                                           element_diameter=2);
    end

    arguments (Repeating)
        u
    end

    arguments
        options.measure_type string = 'dx';
        options.weight = 1;
        options.order = [];
        options.faces = [];     % facet subset for 'ds' (empty = all faces)
    end

    if isempty(options.order)
        % the rule must be exact for the product of all factors: take the
        % max order over every factor (FEFields report their own order;
        % composite nodes report through their signature)
        orders = zeros(1, numel(u));
        for i = 1:numel(u)
            if isa(u{i}, "FEField")
                orders(i) = u{i}.order;
            else
                [~, ~, oi] = read_signature(u{i});
                orders(i) = oi + 1;
            end
        end
        order = max(orders);
    else
        order = options.order;
    end

    u_ = cell(size(u));
    switch options.measure_type
        case 'dx'
            mu = mesh.volume(:);
            [x_local, w] = get_2d_reference_quadrature(order, mesh.ElementType);
            for i = 1:length(u)
                u_{i} = evaluate(u{i}, x_local, mesh); 
            end
        case 'dS'
            mu = mesh.area(mesh.boundaries.');
            [x_local, w] = get_1d_reference_quadrature(order, mesh.ElementType);
            for i = 1:length(u)
                u_{i} = evaluate_boundary(u{i}, x_local, mesh); 
            end
        case 'ds'
            faces = options.faces;
            if isempty(faces)
                faces = 1:size(mesh.faces, 2);
            end
            mu = reshape(mesh.area(faces), [], 1);
            [x_local, w] = get_1d_reference_quadrature(order, mesh.ElementType);
            for i = 1:length(u)
                u__ = evaluate_boundary(u{i}, x_local, mesh);
                u_{i} = element_boundary_to_face(mesh, u__, faces);
            end
        case 'de'
            mu = mesh.area(:);
            [x_local, w] = get_1d_reference_quadrature(order, mesh.ElementType);
            for i = 1:length(u)
                u_{i} = evaluate(u{i}, x_local, mesh); 
            end
    end
    varargout = u_;

    try 
        mu = mu .* options.weight;
    catch
        error('Incorrect weight dimensions.')
    end
end

function u_ = element_boundary_to_face(mesh, u__, faces)
    % gather element-boundary evaluations onto (a subset of) faces using
    % the precomputed incidence table; the table duplicates boundary-face
    % incidences, so both sides are always well defined
    arguments
        mesh AbstractGrid
        u__
        faces = 1:size(mesh.faces, 2);
    end
    K = permute(mesh.incidence(:, 1, faces), [1 3 2]);   % 2 x n_f
    s = permute(mesh.incidence(:, 2, faces), [1 3 2]);
    sz = size(u__, 1:6);
    if sz(5) > 1
        page = reshape(u__, [], sz(5)*sz(6));
        idx = K + sz(5)*(s - 1);                         % (element, side) page
        u_ = reshape(page(:, [idx(1, :), idx(2, :)]), [sz(1:4), numel(faces), 2]);
    else
        page = reshape(u__, [], sz(6));                  % element-independent
        u_ = reshape(page(:, [s(1, :), s(2, :)]), [sz(1:4), numel(faces), 2]);
    end
end