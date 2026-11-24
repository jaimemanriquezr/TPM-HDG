classdef AffineOperator
    % Affine operator of a variational problem: one flat list of terms,
    % each a symbolic IntegralForm array or a pre-assembled IndexedMatrix
    % (e.g. a *Matrices factory output). The condensed/retained split of
    % equations and unknowns is derived at assembly from the spaces'
    % FESpace.Condense, by partitioning the blocks on their [trial, test] ids.
    properties
        Terms cell = {}
    end

    methods
        function obj = AffineOperator(varargin)
            obj.Terms = varargin;
        end

        function obj = plus(obj, term)
            % Operator + term (or term + Operator) appends to the list
            if ~isa(obj, 'AffineOperator')
                [obj, term] = deal(term, obj);
            end
            if isempty(obj)
                obj = AffineOperator();
            end
            obj.Terms{end+1} = term;
        end

        function B = assemble(op, problem)
            % Assemble all terms, merge into one block matrix, and
            % partition into the condensed 2x2 structure:
            %   B.A11 (condensed x condensed)   B.A12 (retained trial)
            %   B.A21 (condensed trial)         B.A22 (retained x retained)
            %   B.A1f / B.A2f  linear blocks per side (0 when absent)
            mesh = problem.Mesh;
            total = [];
            for i = 1:numel(op.Terms)
                t = op.Terms{i};
                if isa(t, 'IntegralForm')
                    all_linear = ~any(arrayfun(@(F) isa(F.integrand, 'BilinearForm'), t));
                    [A, b] = assemble(t, mesh);
                    if ~isa(A, 'IndexedMatrix'); A = IndexedMatrix(A); end
                    total = accumulate(total, A);
                    % An all-linear form returns A == b; skip the duplicate.
                    % Wrap before the double test: IndexedArray subclasses double.
                    if ~all_linear
                        if isa(b, 'IndexedArray'); b = IndexedMatrix(b); end
                        if ~isa(b, 'double')   % bare 0 = no linear terms
                            total = accumulate(total, b);
                        end
                    end
                elseif isa(t, 'IndexedArray')
                    total = accumulate(total, IndexedMatrix(t));
                elseif isa(t, 'IndexedMatrix')
                    total = accumulate(total, t);
                else
                    error('AffineOperator:term', ...
                        'Terms must be IntegralForm arrays or Indexed matrices, got %s.', class(t));
                end
            end

            C = problem.Id.field;    % condensed unknown / equation ids
            R = problem.Id.hybrid;   % retained ids
            trial_id = cellfun(@(c) c(1), total.id);
            test_id  = cellfun(@(c) c(2), total.id);

            % padding data: complete the (trial, test) block layouts per side
            space_dims = arrayfun(@(s) s.Dimension, problem.Space);
            n_el = size(mesh.elements, 2);
            % per-trial face-column count: the sides per element for
            % face-supported retained trials (per-face column structure in
            % dim 7, folded by the local back-substitution), 1 otherwise
            n_sides = size(mesh.boundaries, 1);
            n_face = ones(1, numel(problem.Space));
            n_face(R) = 1 + (n_sides - 1)*([problem.GlobalSpaces.Support] == "face");

            tests_C = unique(test_id(ismember(test_id, C)));
            tests_R = unique(test_id(ismember(test_id, R)));

            B.A11 = pad_ids(total{ismember(trial_id, C) & ismember(test_id, tests_C)}, ...
                            C, tests_C, space_dims, n_el, n_face);
            B.A12 = pad_ids(total{ismember(trial_id, R) & ismember(test_id, tests_C)}, ...
                            R, tests_C, space_dims, n_el, n_face);
            B.A21 = pad_ids(total{ismember(trial_id, C) & ismember(test_id, tests_R)}, ...
                            C, tests_R, space_dims, n_el, n_face);
            B.A22 = pad_ids(total{ismember(trial_id, R) & ismember(test_id, tests_R)}, ...
                            R, tests_R, space_dims, n_el, n_face);
            B.A1f = linear_side(total, trial_id, test_id, C);
            B.A2f = linear_side(total, trial_id, test_id, R);
        end
    end
end

function total = accumulate(total, X)
if isempty(total)
    total = X;
else
    total = total + X;
end
end

function b = linear_side(total, trial_id, test_id, side)
% linear blocks (trial id -1) tested on the given side; a bare 0 when
% absent (downstream sizes missing RHS blocks itself)
mask = trial_id == -1 & ismember(test_id, side);
if any(mask)
    b = total{mask};
else
    b = 0;
end
end

function X = pad_ids(X, trial_ids, test_ids, space_dims, n_el, n_face)
% complete the extracted subset with zero blocks so row/column layouts
% are mutually consistent across the 2x2 partition
present = cell2mat(X.id);
for t = trial_ids(:).'
    for s = test_ids(:).'
        if isempty(present) || ~any(present(:, 1) == t & present(:, 2) == s)
            Z = IndexedArray(zeros(space_dims(s), space_dims(t), 1, 1, n_el, 1, n_face(t)), [t s]);
            X = X + IndexedMatrix(Z);
        end
    end
end
end
