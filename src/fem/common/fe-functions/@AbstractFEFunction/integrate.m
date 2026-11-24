function [I, u_, v_] = integrate(u, v, mesh, options)
%INTEGRATE  Integrate the product of two FE expressions over the mesh.
%
%   [I, u_, v_] = integrate(u, v, mesh, Name = Value, ...)
%
%   u    AbstractFEFunction, or data convertible to Datum.
%   v    AbstractFEFunction (the test side).
%   Options measure_type ('dx' elements, 'dS' faces), order, weight and
%   faces are forwarded to mesh.evaluate_quadrature.
%
%   I  the per-element (per-face for 'dS') integrals as an IndexedArray
%   whose id is [id_u, id_v], read from the signatures; u_, v_ are the
%   evaluated operands.
    arguments
        u
        v AbstractFEFunction
        mesh AbstractGrid = SimplicialGrid(reference_element(), element_diameter=2);
        options.measure_type string = 'dx';
        options.order = [];
        options.weight = 1;
        options.faces = [];
    end
    if ~isa(u, 'AbstractFEFunction')
        try 
            u = Datum(u);
        catch
            error('Invalid first input.')
        end
    end
    op_args = namedargs2cell(options);
    
    switch form_type(u, v, options.measure_type)
        otherwise
            [w, el_measure, u_, v_] = mesh.evaluate_quadrature(u, v, op_args{:});
            el_measure = reshape(el_measure, [1 1 1 1 size(el_measure)]);
            I_ = pagemtimes(w .* v_, 'transpose', u_, 'none');
            if isequal(size(u_, 3:4), size(v_, 3:4))
                % dot product
                I_ = sum(I_, 3:4);
            end
            I = el_measure .* I_;
    end

    if isequal(options.measure_type, 'dS')
        [~, u_dim] = read_signature(u);
        [~, v_dim] = read_signature(v);
        if ~any([u_dim, v_dim] == 1)
            I = sum(I, 6);
        elseif u_dim == 1 && v_dim ~= 1
            I = permute(I, [1 2 3 4 5 7 6]);
        elseif u_dim ~= 1 && v_dim == 1
            % I = permute(I, [1 2 3 4 5 7 6]);
        else
            n_sides = size(mesh.boundaries, 1);
            I_ = zeros([size(I, 1:5) n_sides n_sides]);
            for i = 1:n_sides
                I_(:, :, :, :, :, i, i) = I(:, :, :, :, :, i);
            end
            I = I_;
        end
    end
    [~, ~, ~, id_u] = read_signature(u); 
    [~, ~, ~, id_v] = read_signature(v); 
    I = IndexedArray(I, [id_u, id_v]);
end

function ft = form_type(varargin)
    ft = "None";
end