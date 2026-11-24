% test_interface_fsi.m
% FSI (PEERS stress-velocity fluid over PEERS linear elasticity) on
% InterfaceProblem, compared dof-for-dof (tol 1e-10) against
% FluidStructureCouplingOracle, the frozen reference implementation, for all
% schemes: steady, pinned, onestep at time 0, onestep marching and Newmark,
% with nonzero volume and skeleton histories. Polynomial manufactured data,
% so both volume-source rules integrate exactly. Also checks the symbolic
% transmission rows against fsi_rows under the pinned closure.
clc; clear;

y0 = 0.5;
interface = FlatCurveHorizontal(y0);
mu_f = 1.0;
E_s = 1.0;
nu_s = 0.25;
dt = 1e-2;   % mass/dt^2 dominates the local systems; 1e-2 keeps their
             % conditioning ~1e4 so roundoff stays below the 1e-10 bound
rho_f = 1.0;
rho_s = 1.0;
beta_nm = 1/4;
gamma_nm = 1/2;

n_fail = 0;
for scheme = ["steady", "pinned", "onestep", "marching", "newmark"]
    % "pinned" = steady under the pinned pressure closure (no ConstantBasis
    % multiplier; the alpha scalar is strong-pinned, oracle's p0_h = 0 path)
    oracle_scheme = scheme;
    if scheme == "marching", oracle_scheme = "onestep"; end
    if scheme == "pinned",   oracle_scheme = "steady";  end
    closure = "constrained";
    if scheme == "pinned", closure = "pinned"; end
    for k = [1, 2]
        d1 = k + 1;
        d2 = nchoosek(k + 2, 2);
        Pk = DubinerBasis(k);

        [fluid_pb2, solid_pb2, interface_sources, problem_data] = ...
            FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
            "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
            0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = closure);

        for h = [1/8, 1/16]
            d_gap = h^2;
            [fluid_mesh, solid_mesh] = fitted_split(interface, [h, h], d_gap);
            nel_F = size(fluid_mesh.elements, 2);
            nel_S = size(solid_mesh.elements, 2);
            nF_S  = size(solid_mesh.faces, 2);

            % deterministic nonzero histories (any data must match)
            wave = @(n, m, c) reshape(sin(c + (1:n*m)), n, m);
            u0_F  = wave(2*d2, nel_F, 0.1);   % fluid initial velocity (volume)
            a0_F  = wave(2*d2, nel_F, 0.2);   % fluid u_t history (theta < 1)
            u0_S  = wave(2*d2, nel_S, 0.3);   % solid initial displacement
            v0_S  = wave(2*d2, nel_S, 0.4);
            a0_S  = wave(2*d2, nel_S, 0.5);
            u00_S = wave(2*d2, nel_S, 0.6);   % previous displacement (marching)
            u0h_S = wave(2*d1, nF_S, 0.7);    % skeleton histories (T2 row)
            v0h_S = wave(2*d1, nF_S, 0.8);
            a0h_S = wave(2*d1, nF_S, 0.9);

            % per-scheme configuration: the state (factory options and the
            % oracle's Data.state), the T2 history of InterfaceProblem.fsi
            % and the oracle's skeleton Data fields
            fluid_opts = {};  solid_opts = {};  coupling_opts = {};  oracle_opts = {};  rows_opts = {};
            fluid_pb2.Data.density = 0.0;
            solid_pb2.Data.density = 0.0;
            fluid_pb2.Data.state = struct.empty;
            solid_pb2.Data.state = struct.empty;
            switch scheme
                case {"steady", "pinned"}
                    % densities stay 0; cdt = dt, no history
                case "onestep"
                    fluid_pb2.Data.density = rho_f;
                    solid_pb2.Data.density = rho_s;
                    solid_pb2.Data.time = 0.0;
                    state = struct("fluid", struct("velocity", u0_F), ...
                        "solid", struct("displacement", u0_S, "velocity", v0_S));
                    history = struct("displacement_skeleton", u0h_S);
                    fluid_opts = {'density', rho_f, 'time_step', dt};
                    solid_opts = {'density', rho_s, 'time_step', dt};
                case "marching"
                    fluid_pb2.Data.density = rho_f;
                    solid_pb2.Data.density = rho_s;
                    solid_pb2.Data.time = dt;
                    state = struct("fluid", struct("velocity", u0_F), ...
                        "solid", struct("displacement", u0_S, "displacement_previous", u00_S));
                    history = struct("displacement_skeleton", u0h_S);
                    fluid_opts = {'density', rho_f, 'time_step', dt};
                    solid_opts = {'density', rho_s, 'time_step', dt};
                case "newmark"
                    fluid_pb2.Data.density = rho_f;
                    fluid_pb2.Data.theta = 0.5;      % Crank-Nicolson
                    solid_pb2.Data.density = rho_s;
                    solid_pb2.Data.time_scheme = "newmark";
                    solid_pb2.Data.newmark_beta = beta_nm;
                    solid_pb2.Data.newmark_gamma = gamma_nm;
                    state = struct("fluid", struct("velocity", u0_F, "acceleration", a0_F), ...
                        "solid", struct("displacement", u0_S, "velocity", v0_S, "acceleration", a0_S));
                    history = struct("displacement_skeleton", u0h_S, ...
                        "velocity_skeleton", v0h_S, "acceleration_skeleton", a0h_S);
                    fluid_opts = {'density', rho_f, 'time_step', dt, 'theta', 0.5};
                    solid_opts = {'density', rho_s, 'time_step', dt, ...
                        'time_scheme', "newmark", 'newmark_beta', beta_nm};
                    % a nonzero history extension too: it is per patch, so
                    % its width comes from the induced partition
                    [~, part_S] = get_partitions(interface, fluid_mesh, solid_mesh);
                    hx = wave(2*d1, size(part_S.faces, 3), 1.1);
                    % on the oracle and the fsi_rows path only: the
                    % symbolic transmission below has no history transport
                    oracle_opts = {'history_extension', hx};
                    coupling_opts = {'newmark_beta', beta_nm, 'newmark_gamma', gamma_nm};
            end
            if ~ismember(scheme, ["steady", "pinned"])
                fluid_opts = [fluid_opts, {'state', state.fluid}]; %#ok<AGROW>
                solid_opts = [solid_opts, {'state', state.solid}]; %#ok<AGROW>
                % the oracle reads the volume histories through Data.state
                % and the skeleton histories from Data fields
                fluid_pb2.Data.state = state;
                solid_pb2.Data.state = state;
                for f = string(fieldnames(history)).'
                    solid_pb2.Data.("initial_" + f) = history.(f);
                end
                coupling_opts = [coupling_opts, {'history', history}]; %#ok<AGROW>
                if scheme == "newmark"
                    % the later history option wins: the same one plus hx
                    history.extension = hx;
                    rows_opts = {'history', history};
                end
            end

            %% frozen reference implementation
            [fs2, ss2] = FluidStructureCouplingOracle.solve(fluid_pb2, solid_pb2, ...
                fluid_mesh, solid_mesh, interface, ...
                interface_sources{1}, interface_sources{2}, oracle_scheme, ...
                'Order', k, 'quiet', true, oracle_opts{:});

            %% InterfaceProblem
            fluid_gspaces = [LegendreBasis(k).to_vector(2), DubinerBasis(0)];
            if closure == "constrained"
                fluid_gspaces = [fluid_gspaces, ConstantBasis()];
            end
            fluid_pb1 = FEProblem( ...
                Spaces  = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), FESpace(fluid_gspaces)], ...
                Order = k);
            fluid_pb1.Mesh = fluid_mesh;
            [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, fluid_mesh, mu_f, ...
                problem_data.f_F, "pressure_mean", closure, fluid_opts{:});
            fluid_pb1.Operator = AffineOperator(lhs + rhs);
            fluid_pb1.Operator = fluid_pb1.Operator + (flux);
            fluid_pb1.BoundaryConditions = StrongEssentialBoundaryConditions( ...
                Datum(problem_data.u_F), "Dirichlet", projection_order = k);

            solid_pb1 = FEProblem( ...
                Spaces  = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense=true), FESpace(LegendreBasis(k).to_vector(2))], ...
                Order = k);
            solid_pb1.Mesh = solid_mesh;
            [solid_lhs, solid_flux, solid_rhs] = ...
                LinearElasticityMatrices(k, solid_mesh, E_s, nu_s, solid_opts{:});
            solid_pb1.Operator = AffineOperator(solid_lhs);
            solid_pb1.Operator = solid_pb1.Operator + (solid_flux);
            if ~isempty(solid_rhs)
                solid_pb1.Operator = solid_pb1.Operator + solid_rhs;
            end
            w = TestFunction(Pk);
            F1 = double(Datum(problem_data.f_S{1}) * w * Measure(solid_mesh, 'dx'));
            F2 = double(Datum(problem_data.f_S{2}) * w * Measure(solid_mesh, 'dx'));
            F = IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));
            solid_pb1.Operator = solid_pb1.Operator + F;
            solid_pb1.BoundaryConditions = StrongEssentialBoundaryConditions( ...
                Datum(problem_data.u_S), "Dirichlet", projection_order = k);

            coupled = InterfaceProblem.fsi(fluid_pb1, solid_pb1, interface, ...
                'viscosity', mu_f, 'young_modulus', E_s, 'poisson_ratio', nu_s, ...
                'domain_measure', problem_data.measure_F, ...
                'scheme', oracle_scheme, 'time_step', dt, ...
                'source_neumann', interface_sources{1}, ...
                'source_dirichlet', interface_sources{2}, ...
                'pressure_mean', closure, coupling_opts{:}, rows_opts{:});
            [fs1, ss1, gs] = solve(coupled, k);

            %% symbolic transmission, pinned closure: the coupled traction
            %% balance fixes the fluid pressure level
            if scheme == "pinned" || scheme == "newmark"
                fp = FEProblem( ...
                    Spaces  = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])], ...
                    Order = k);
                fp.Mesh = fluid_mesh;
                [lhsP, fluxP, rhsP] = IncompressibleElasticityMatrices(k, fluid_mesh, mu_f, ...
                    problem_data.f_F, "pressure_mean", "pinned", fluid_opts{:});
                fp.Operator = AffineOperator(lhsP + rhsP);
                fp.Operator = fp.Operator + (fluxP);
                fp.BoundaryConditions = StrongEssentialBoundaryConditions( ...
                    Datum(problem_data.u_F), "Dirichlet", projection_order = k);
                sp = solid_pb1;
                if scheme == "pinned"
                    fsP = fs1; ssP = ss1;   % the rows path is already pinned
                else
                    coupledP = InterfaceProblem.fsi(fp, sp, interface, ...
                        'viscosity', mu_f, 'young_modulus', E_s, 'poisson_ratio', nu_s, ...
                        'domain_measure', problem_data.measure_F, ...
                        'scheme', oracle_scheme, 'time_step', dt, ...
                        'source_neumann', interface_sources{1}, ...
                        'source_dirichlet', interface_sources{2}, ...
                        'pressure_mean', "pinned", coupling_opts{:});
                    [fsP, ssP] = solve(coupledP, k);
                end
                [sg1, u1v, ~, ~, uhat1, ~] = get_trial_functions(fp);
                [~, ~, ~, ~, mu1v, ~] = get_test_functions(fp);
                [sg2, u2v, ~, uhat2] = get_trial_functions(sp);
                [~, ~, ~, mu2v] = get_test_functions(sp);
                gm1 = get_nth_trial(fp, 3); gm2 = get_nth_trial(sp, 3);
                n1_ = InterfaceNormal(fluid_mesh); n2_ = InterfaceNormal(solid_mesh);
                tau_ = Stabilization();
                dI1 = Measure(fluid_mesh, "dI"); dI2 = Measure(solid_mesh, "dI");
                c11 = 1/(2*mu_f); c21 = 1/(4*mu_f);
                [lam_s, mus_] = LinearElasticityProblem.get_lame_parameters(E_s, nu_s);
                c12 = 1/(2*mus_); c22 = lam_s/(2*mus_*(2*lam_s + 2*mus_));
                Jm = [0, 1; -1, 0];
                T1_ = Datum(interface_sources{1}); T2_ = Datum(interface_sources{2});
                if scheme == "pinned"
                    cdt = dt; h0v = [];
                else
                    cdt = (beta_nm/gamma_nm)*dt;
                    h0v = u0h_S + dt*(1 - beta_nm/gamma_nm)*v0h_S ...
                        + dt^2*((1/2 - beta_nm) - beta_nm*(1 - gamma_nm)/gamma_nm)*a0h_S;
                end
                T1s = ( Extrapolate(sg1)*n1_ - (mu_f*tau_)*(u1v - uhat1) ...
                      + Extrapolate(sg2)*n2_ - tau_*(u2v - uhat2) ...
                      - T1_ ) * mu1v * dI1;
                T2s = cdt*uhat1 ...
                    + cdt*PathIntegral(c11*sg1 - c21*trace(sg1)*eye(2) + gm1*Jm) ...
                    - uhat2 - PathIntegral(c12*sg2 - c22*trace(sg2)*eye(2) + gm2*Jm) ...
                    - cdt*T2_;
                if ~isempty(h0v)
                    H0_ = FEField(FESpace(LegendreBasis(k).to_vector(2), solid_mesh), h0v);
                    H0_.mesh = solid_mesh;
                    T2s = T2s + H0_;
                end
                T2s = T2s * mu2v * dI2;
                coupledS = InterfaceProblem([fp, sp], interface, struct("kind", "none"));
                coupledS.transmission = [T1s, T2s];
                [fsS, ssS] = solve(coupledS, k);
                d_symb = 0;
                for fld = 1:3
                    d_symb = max(d_symb, max(abs(fsS.Unknowns.field(fld).dof ...
                                              - fsP.Unknowns.field(fld).dof), [], 'all'));
                    d_symb = max(d_symb, max(abs(ssS.Unknowns.field(fld).dof ...
                                              - ssP.Unknowns.field(fld).dof), [], 'all'));
                end
                d_symb = max(d_symb, max(abs(fsS.Unknowns.hybrid(1).dof ...
                                          - fsP.Unknowns.hybrid(1).dof), [], 'all'));
                d_symb = max(d_symb, max(abs(ssS.Unknowns.hybrid(1).dof ...
                                          - ssP.Unknowns.hybrid(1).dof), [], 'all'));
                fprintf("%-8s k = %d, h = %.4f: symbolic vs fsi_rows (pinned) max dof diff = %.2e\n", ...
                    scheme, k, h, d_symb);
                if d_symb > 1e-10
                    fprintf("%s k = %d symbolic transmission FAILED\n", scheme, k);
                    n_fail = n_fail + 1;
                end
            end

            %% dof-for-dof comparison
            dsF = max(abs(fs1.Unknowns.field(1).dof - fs2.Unknowns.field(1).dof), [], "all");
            duF = max(abs(fs1.Unknowns.field(2).dof - fs2.Unknowns.field(2).dof), [], "all");
            dgF = max(abs(fs1.Unknowns.field(3).dof - fs2.Unknowns.field(3).dof), [], "all");
            uhatF1 = reshape(fs1.Unknowns.hybrid(1).dof, 2*d1, []);
            duhF = max(abs(uhatF1 - fs2.Unknowns.hybrid(1).dof), [], "all");
            drho = max(abs(reshape(fs1.Unknowns.hybrid(2).dof, 1, []) ...
                         - reshape(fs2.Unknowns.hybrid(2).dof, 1, [])), [], 'all');

            dsS = max(abs(ss1.Unknowns.field(1).dof - ss2.Unknowns.field(1).dof), [], "all");
            duS = max(abs(ss1.Unknowns.field(2).dof - ss2.Unknowns.field(2).dof), [], "all");
            dgS = max(abs(ss1.Unknowns.field(3).dof - ss2.Unknowns.field(3).dof), [], "all");
            uhatS1 = reshape(ss1.Unknowns.hybrid(1).dof, 2*d1, []);
            duhS = max(abs(uhatS1 - ss2.Unknowns.hybrid(1).dof), [], "all");

            worst = max([dsF, duF, dgF, duhF, drho, dsS, duS, dgS, duhS]);
            fprintf("%-8s k = %d, h = %.4f: worst dof diff = %.3e (sF %.1e, uF %.1e, gF %.1e, uhatF %.1e, rho %.1e, sS %.1e, uS %.1e, gS %.1e, uhatS %.1e)\n", ...
                scheme, k, h, worst, dsF, duF, dgF, duhF, drho, dsS, duS, dgS, duhS);
            if worst > 1e-10
                fprintf("%s k = %d, h = %.4f ORACLE COMPARISON FAILED\n", scheme, k, h);
                n_fail = n_fail + 1;
            end
        end
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_fsi:failed", "%d checks failed", n_fail);
end

function t = get_nth_trial(pb, n)
out = cell(1, numel(pb.Space));
[out{:}] = get_trial_functions(pb);
t = out{n};
end
