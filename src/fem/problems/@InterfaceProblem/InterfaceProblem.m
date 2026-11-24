classdef InterfaceProblem
    % TPM interface problem over two FEProblems meeting at an AbstractCurve.
    %
    % Each subproblem carries its own mesh and global-support layout. The
    % coupled system concatenates the two condensed systems
    % (condense_global_system), deletes each side's interface-face flux
    % rows and substitutes the transmission rows of the coupling.
    %
    % Coupling factories are static methods (one per pairing):
    %   InterfaceProblem.darcy_darcy(problem_1, problem_2, interface, ...)
    %   InterfaceProblem.stokes_stokes(problem_1, problem_2, interface, ...)
    %   InterfaceProblem.stokes_biot(fluid_problem, biot_problem, interface, ...)
    properties
        problems  (1, 2) FEProblem
        % one AbstractCurve, or a per-side pair {curve_1, curve_2} on one
        % parameter (interface_curves; fsi coupling only), validated by
        % the constructor
        interface
        coupling  struct   % kind + physical parameters, set by the factories
        options   struct   % quadrature orders, vertex tolerances, ...
        % symbolic transmission rows (ExtendedFunction array, each with a
        % test and a Measure(mesh_i, "dI")); when non-empty, solve
        % assembles these instead of the built-in *_rows
        transmission ExtendedFunction = ExtendedFunction.empty
    end

    methods
        function obj = InterfaceProblem(problems, interface, coupling, options)
            arguments
                problems  (1, 2) FEProblem
                interface {mustBeInterfaceCurves}
                coupling  struct = struct("kind", "none")
                options.quadrature_patch double = []   % default k+1
                options.quadrature_path  double = []   % default k+2
                options.vertex_tol       (1,1) double = 1e-8
                options.vertex_tol_relative double = []
                % faces with a zero-length image under the closest-point
                % map (a sharp tip of Sigma): see get_partitions. "collapse"
                % is a fallback where "spread" fails (orphanFace); where
                % "spread" works it is the more accurate choice
                options.tip_mapping (1,1) string ...
                    {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
                options.collapsed_normal (1,1) string ...
                    {mustBeMember(options.collapsed_normal, ["path", "face", "bisector"])} = "path"
                % per-side point-to-interface maps (get_partitions); ""
                % keeps the curve's own
                options.mapping_1 (1,1) string = ""
                options.mapping_2 (1,1) string = ""
            end
            obj.problems  = problems;
            obj.interface = interface;
            obj.coupling  = coupling;
            obj.options   = options;
        end

        [solution_1, solution_2, global_system, cache] = solve(obj, k, options)
    end

    methods (Static)
        function obj = darcy_darcy(problem_1, problem_2, interface, kappa_1, kappa_2, options)
            % Darcy-Darcy transmission on Sigma:
            %   [[p]] = 0 imposed on the transferred traces (T2 row),
            %   [[q . n]] = 0 on the extrapolated numerical fluxes (T1 row).
            % Each problem must carry local_spaces = [VPk, Pk] and a scalar
            % trace global space (the DarcyMatrices layout).
            arguments
                problem_1 (1,1) FEProblem
                problem_2 (1,1) FEProblem
                interface AbstractCurve
                kappa_1 (2,2) double = eye(2)
                kappa_2 (2,2) double = eye(2)
                options.quadrature_patch double = []
                options.quadrature_path  double = []
                options.vertex_tol       (1,1) double = 1e-8
                options.vertex_tol_relative double = []
            end
            coupling = struct("kind", "darcy_darcy", ...
                "kappa", {{kappa_1, kappa_2}});
            opts = namedargs2cell(options);
            obj = InterfaceProblem([problem_1, problem_2], interface, coupling, opts{:});
        end

        function obj = stokes_stokes(problem_1, problem_2, interface, nu_1, nu_2, options)
            % Stokes-Stokes transmission between two velocity gradient-
            % velocity-pressure (L-u-p) subdomains, set as symbolic
            % transmission forms (obj.transmission):
            %   T1 (tested on side-1 interface faces): traction balance on
            %       the extrapolated fields, stabilization pair on the face,
            %       sum_i < nu_i E(L_i) n_i - E(p_i) n_i
            %               - nu_i tau (u_i - u_hat_i), mu_1 > = 0
            %   T2 (tested on side-2 interface faces): transferred
            %       velocity continuity,
            %       < (u_hat_1 + int_path L_1 . m) - (u_hat_2 + int_path L_2 . m), mu_2 > = 0
            % with mu_i the side-i trace test functions and m the transfer
            % path vector. Each problem must carry its mesh and the
            % StokesMatrices layout with kinematic viscosity nu_i: local
            % [L, u, p, xi], global [u_hat, rho, mu_c] (mu_c = the
            % ConstantBasis zero-mean multiplier). At least one side must
            % carry this pressure closure; the other may omit mu_c (and
            % its two blocks), its pressure level then being fixed through
            % the traction balance T1.
            arguments
                problem_1 (1,1) FEProblem
                problem_2 (1,1) FEProblem
                interface AbstractCurve
                nu_1 (1,1) double = 1.0
                nu_2 (1,1) double = 1.0
                options.quadrature_patch double = []
                options.quadrature_path  double = []
                options.vertex_tol       (1,1) double = 1e-8
                options.vertex_tol_relative double = []
            end
            assert(~isempty(problem_1.Mesh) && ~isempty(problem_2.Mesh), ...
                "InterfaceProblem:mesh", "both subproblems must carry their meshes");
            opts = namedargs2cell(options);
            obj = InterfaceProblem([problem_1, problem_2], interface, opts{:});

            [L1, u1, p1, ~, uhat1] = get_trial_functions(problem_1);
            [L2, u2, p2, ~, uhat2] = get_trial_functions(problem_2);
            [~, ~, ~, ~, mu1] = get_test_functions(problem_1);
            [~, ~, ~, ~, mu2] = get_test_functions(problem_2);
            n1_ = InterfaceNormal(problem_1.Mesh);
            n2_ = InterfaceNormal(problem_2.Mesh);
            tau_ = Stabilization();
            dI1 = Measure(problem_1.Mesh, "dI");
            dI2 = Measure(problem_2.Mesh, "dI");

            T1 = ( nu_1*(Extrapolate(L1)*n1_) - Extrapolate(p1)*n1_ ...
                     - (nu_1*tau_)*(u1 - uhat1) ...
                 + nu_2*(Extrapolate(L2)*n2_) - Extrapolate(p2)*n2_ ...
                     - (nu_2*tau_)*(u2 - uhat2) ) * mu1 * dI1;
            T2 = ( uhat1 + PathIntegral(L1) ...
                 - uhat2 - PathIntegral(L2) ) * mu2 * dI2;
            obj.transmission = [T1, T2];
        end

        function obj = stokes_darcy(stokes_problem, darcy_problem, interface, options)
            % Stokes-Darcy transmission with Beavers-Joseph-Saffman
            % friction: mass conservation, BJS and pressure balance, plus
            % the extrapolated-pressure scalar alpha as a coupling-owned
            % unknown. Side 1 = velocity-gradient Stokes (StokesMatrices
            % layout, constrained closure), side 2 = Darcy (DarcyMatrices
            % layout).
            arguments
                stokes_problem (1,1) FEProblem
                darcy_problem  (1,1) FEProblem
                interface AbstractCurve
                options.viscosity      (1,1) double = 1.0
                options.diffusion      (2,2) double = eye(2)
                options.friction_bjs   (1,1) double = 0.0
                options.domain_measure (1,1) double = 1.0
                options.source_neumann   cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
                options.source_dirichlet (1,1) function_handle = @(x,y) 0*x.*y
                options.quadrature_patch double = []
                options.quadrature_path  double = []
                options.vertex_tol       (1,1) double = 1e-8
                options.vertex_tol_relative double = []
            end
            coupling = struct("kind", "stokes_darcy", "n_extra", 1, ...
                "viscosity", options.viscosity, ...
                "diffusion", options.diffusion, ...
                "friction_bjs", options.friction_bjs, ...
                "domain_measure", options.domain_measure, ...
                "source_neumann", {options.source_neumann}, ...
                "source_dirichlet", options.source_dirichlet);
            obj = InterfaceProblem([stokes_problem, darcy_problem], interface, ...
                coupling, ...
                quadrature_patch = options.quadrature_patch, ...
                quadrature_path = options.quadrature_path, ...
                vertex_tol = options.vertex_tol, ...
                vertex_tol_relative = options.vertex_tol_relative);
        end

        function obj = fsi(fluid_problem, solid_problem, interface, options)
            % Fluid-structure transmission: side 1 = stress-velocity PEERS
            % fluid (IncompressibleElasticityMatrices layout), side 2 =
            % PEERS linear elasticity (LinearElasticityMatrices layout).
            %   T1: traction balance on the extrapolated stresses (both
            %       sides), alpha = the fluid pressure-level scalar as a
            %       coupling-owned extra unknown (sliver-weighted mean row);
            %   T2: kinematic row < cdt*u_F - u_S, mu_2 > on the
            %       transferred velocities/displacements. Under
            %       interface_law = "navier"/"full" only its normal
            %       projection is imposed and the tangential rows carry
            %       the fluid Navier law instead.
            %       cdt and the history h0 encode the scheme:
            %       steady  cdt = dt, h0 = 0 (one backward-Euler step
            %               from rest: < dt*u_F - u_S, mu_2 >);
            %       steady, time_step = Inf: the stationary problem, the
            %               row divided by dt in the limit dt -> Inf,
            %               < u_F, mu_2 > = < T2, mu_2 > (the solid velocity
            %               vanishes; no displacement term, no history;
            %               the pressure level is then fixed by the
            %               "constrained" closure, which is required);
            %       onestep cdt = dt, h0 = u0 (backward Euler);
            %       newmark cdt = (beta/gamma)*dt,
            %               h0 = u0 + dt(1-beta/gamma) v0
            %                  + dt^2[(1/2-beta) - beta(1-gamma)/gamma] a0.
            %       u0, v0, a0 are the solid skeleton histories of the
            %       history option; its transport Ext(h0) onto Sigma is
            %       added to h0 after the per-patch gather (see the
            %       history option).
            % Volume time terms belong to the subproblems' local matrices,
            % not to the coupling. Closed-loop anchors (circulation, rigid
            % modes, traction resultants) are not supported.
            %
            % Per-side curves (Lagrangian solid): interface may be the pair
            % {curve_fluid, curve_solid} on one parameter s, the solid
            % reference Sigma_0, s -> X(s), and the current interface
            % Sigma(t), s -> x(s) (e.g. DeformedCurve: x = X + u_S(X)).
            % The fluid paths end on Sigma(t), the solid paths on Sigma_0,
            % and patches pair through s (get_partitions, interface_2).
            % Each side's normal is its own curve's. Both rows keep the
            % per-patch measure of the one-curve coupling (the tested
            % face's length per unit s, i.e. ds up to a per-face row
            % scaling), at the same s-nodes; the pointwise Nanson factor
            % J(s) = |x'(s)|/|X'(s)| then writes the traction balance on
            % Sigma_0:
            %   T1: J*(sigma_F n_F - mu*tau_F*(u_F - u_hat_F) + p0 n_F)
            %       + sigma_S n_S - tau_S*(u_S - u_hat_S) = T1, i.e. every
            %       fluid term (extrapolated stress, its stabilization
            %       pair and the alpha column) carries J; the solid terms
            %       and the datum do not;
            %   T2: unchanged, u_F at x(s) against u_S at X(s), no J (a
            %       pointwise equality at the same s).
            % source_neumann is evaluated at the fluid points x(s) and
            % source_dirichlet at the solid points X(s). The sliver of
            % the alpha row is the fluid's, against Sigma(t). With one
            % curve J = 1 and nothing is multiplied. The pair supports
            % interface_law = "noslip" only.
            arguments
                fluid_problem (1,1) FEProblem
                solid_problem (1,1) FEProblem
                % one curve, or the pair {curve_fluid, curve_solid}: see
                % "Per-side curves" above
                interface {mustBeInterfaceCurves}
                options.viscosity      (1,1) double = 1.0
                options.young_modulus  (1,1) double = 1.0
                options.poisson_ratio  (1,1) double = 0.25
                options.domain_measure (1,1) double = 1.0
                options.scheme (1,1) string ...
                    {mustBeMember(options.scheme, ["steady", "onestep", "newmark"])} = "steady"
                options.time_step      (1,1) double = 1.0
                options.newmark_beta   (1,1) double = 1/4
                options.newmark_gamma  (1,1) double = 1/2
                % The T2 history (onestep / newmark), a struct with fields
                % named as in the solid state (the solid side of
                % FluidStructureCoupling.initial_state, which may be
                % passed as is); [], 0 or a missing field = zero:
                %   displacement_skeleton, velocity_skeleton,
                %   acceleration_skeleton  the solid skeleton traces u0,
                %       v0, a0 (2*(k+1) x n_faces)
                %   stress, rotation (onestep and newmark),
                %   stress_velocity, rotation_velocity,
                %   stress_acceleration, rotation_acceleration (newmark)
                %       the solid PEERS stress and rotation dofs of u0,
                %       v0, a0, from which solve builds the transport
                %       Ext(h0) (trace_extension_coefficients on the
                %       solid side of its partition; under newmark from
                %       the same combination as h0 below). The kinematic
                %       row compares a path-extended trace at n+1 against
                %       the face trace of h0; on a retreated mesh that is
                %       an O(delta/dt) datum error, and adding Ext(h0)
                %       after the per-patch gather puts both sides on
                %       Sigma. No stress history: no transport.
                %   extension  Ext(h0) given directly, (2*d1) x n_patches
                %       coefficients; when non-empty it is used instead of
                %       the one built from the stress history.
                % struct.empty = zero history, no transport.
                options.history struct = struct.empty
                options.source_neumann   cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
                options.source_dirichlet cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
                options.pressure_mean (1,1) string ...
                    {mustBeMember(options.pressure_mean, ["constrained", "pinned"])} = "constrained"
                options.patch_weight (1,1) string {mustBePatchWeight} = "rectangular"
                options.sliver_correction (1,1) logical = true
                % Interface law. "noslip" is the classical coupling
                % (T1 traction balance + T2 full velocity continuity).
                % "navier" / "full" follow Cesik, Kampschulte &
                % Schwarzacher (arXiv:2603.12030): the kinematic coupling
                % reduces to impermeability, (v - u_S).n = 0, and the
                % tangential rows carry the fluid Navier law
                %   t.(sigma_F n) + a t.(cdt*v - u_S) = t.g,
                % a = slip_coefficient >= 0 ("full" is a = 0). The friction
                % acts on the scheme-weighted velocity difference of the T2
                % row, so the units of a depend on the time scheme.
                % Verified at dt = 1 only; check a when dt != 1.
                % T1 is unchanged; only the fluid traction is prescribed.
                % The no-slip coupling is the a -> infinity limit.
                options.interface_law (1,1) string ...
                    {mustBeMember(options.interface_law, ["noslip", "navier", "full"])} = "noslip"
                % used by "navier" only; nonzero with another law is an error
                options.slip_coefficient (1,1) double {mustBeNonnegative} = 0.0
                % data for the Navier row: its tangential projection is the
                % prescribed fluid tangential traction g = sigma_F n_F on Sigma
                options.source_slip cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
                % Test the traction balance (T1) with the solid's interface
                % traces and the kinematic row (T2) with the fluid's, instead
                % of the other way round; for interface_law "noslip" only.
                % The kinematic row then has as many rows as the fluid has
                % interface traces, which keeps it from over-constraining a
                % coarse fluid against a finer solid.
                options.swap_transmission (1,1) logical = false
                options.quadrature_patch double = []
                options.quadrature_path  double = []
                options.vertex_tol       (1,1) double = 1e-8
                options.vertex_tol_relative double = []
                % see get_partitions; a collapsed face carries the row
                % tested on its own side over its whole length and none
                % of the row tested on the other side. solve requires the
                % kinematic row on the collapsed side (fluid faces:
                % swap_transmission = true; solid faces: false) and
                % interface_law = "noslip"; collapsed_normal then does
                % not enter the system.
                % Accuracy: "collapse" is a fallback that makes meshes
                % solvable where "spread" fails (orphanFace when the tied
                % run straddles the parameter seam, as on an exact cusp).
                % Where "spread" works it is the more accurate choice: at
                % k = 2 on meshes with six tied faces at a cusp, collapse
                % gave errors about 4-5 times larger, and with collapse
                % the error near a cusp depends on the tie pattern and is
                % not observed to converge for k >= 2. Do not use it for
                % convergence studies.
                options.tip_mapping (1,1) string ...
                    {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
                options.collapsed_normal (1,1) string ...
                    {mustBeMember(options.collapsed_normal, ["path", "face", "bisector"])} = "path"
                options.mapping_1 (1,1) string = ""
                options.mapping_2 (1,1) string = ""
            end
            if isinf(options.time_step)
                % Without the solid term the kinematic rows no longer carry
                % the volume constraint (the Sigma-integral of the row) that
                % fixes the fluid pressure level, so the closure must.
                assert(options.scheme == "steady" && options.time_step > 0 ...
                    && options.pressure_mean == "constrained", ...
                    "InterfaceProblem:infiniteTimeStep", ...
                    "time_step = Inf (the stationary problem) requires " + ...
                    "scheme = ""steady"" and pressure_mean = ""constrained"".");
            end
            assert(~options.swap_transmission || options.interface_law == "noslip", ...
                "InterfaceProblem:swapTransmission", ...
                "swap_transmission is implemented for interface_law = ""noslip"" only.");
            % pressure_mean must agree with the fluid factory's closure
            % (constrained <=> the ConstantBasis multiplier is present)
            has_mu_c = any(arrayfun(@(sp) isa(sp.Basis, 'ConstantBasis'), ...
                fluid_problem.GlobalSpaces));
            assert(has_mu_c == (options.pressure_mean == "constrained"), ...
                "InterfaceProblem:pressure_mean", ...
                "pressure_mean = ""%s"" but the fluid problem was built " + ...
                "with the %s closure (ConstantBasis %s in global_spaces).", ...
                options.pressure_mean, ...
                ternary(has_mu_c, "constrained", "pinned"), ...
                ternary(has_mu_c, "present", "absent"));
            dt = options.time_step;
            h0 = 0.0;
            % (sigma, gamma) of h0, from which solve builds Ext(h0);
            % [] = no stress history
            history_stress = [];
            history_rotation = [];
            % weight of the solid trace (and its history) in the kinematic
            % row; 0 only for the stationary limit time_step = Inf
            solid_weight = 1.0;
            switch options.scheme
                case "steady"
                    cdt = dt;
                    if isinf(dt)
                        cdt = 1.0;
                        solid_weight = 0.0;
                    end
                case "onestep"
                    cdt = dt;
                    h0 = history_field(options.history, "displacement_skeleton");
                    history_stress   = history_field(options.history, "stress");
                    history_rotation = history_field(options.history, "rotation");
                case "newmark"
                    beta_nm  = options.newmark_beta;
                    gamma_nm = options.newmark_gamma;
                    cdt = (beta_nm/gamma_nm)*dt;
                    c_v = dt*(1 - beta_nm/gamma_nm);
                    c_a = dt^2*((1/2 - beta_nm) - beta_nm*(1 - gamma_nm)/gamma_nm);
                    h0 = history_field(options.history, "displacement_skeleton") ...
                       + c_v * history_field(options.history, "velocity_skeleton") ...
                       + c_a * history_field(options.history, "acceleration_skeleton");
                    history_stress = history_field(options.history, "stress") ...
                       + c_v*history_field(options.history, "stress_velocity") ...
                       + c_a*history_field(options.history, "stress_acceleration");
                    history_rotation = history_field(options.history, "rotation") ...
                       + c_v*history_field(options.history, "rotation_velocity") ...
                       + c_a*history_field(options.history, "rotation_acceleration");
            end
            if isequal(history_stress, 0.0) && isequal(history_rotation, 0.0)
                % zero stress history: no transport
                history_stress = [];
                history_rotation = [];
            end
            history_extension = [];
            if isfield(options.history, "extension")
                history_extension = options.history.extension;
            end
            if ~isempty(history_extension)
                % the explicit transport is used instead
                history_stress = [];
                history_rotation = [];
            end
            if options.interface_law ~= "navier" && options.slip_coefficient ~= 0
                error("InterfaceProblem:slipCoefficientUnused", ...
                    "slip_coefficient = %g is ignored by interface_law = ""%s""; " + ...
                    "use ""navier"" for a > 0 (""full"" IS a = 0).", ...
                    options.slip_coefficient, options.interface_law);
            end
            if options.interface_law == "navier" && options.slip_coefficient == 0
                error("InterfaceProblem:slipCoefficientZero", ...
                    "interface_law = ""navier"" with slip_coefficient = 0 is " + ...
                    "full slip; ask for ""full"" so the intent is explicit.");
            end
            coupling = struct("kind", "fsi", "n_extra", 1, ...
                "pin_extra", options.pressure_mean == "pinned", ...
                "viscosity", options.viscosity, ...
                "young_modulus", options.young_modulus, ...
                "poisson_ratio", options.poisson_ratio, ...
                "domain_measure", options.domain_measure, ...
                "cdt", cdt, "solid_weight", solid_weight, ...
                "history_skeleton", h0, ...
                "history_extension", history_extension, ...
                "history_stress", history_stress, ...
                "history_rotation", history_rotation, ...
                "patch_weight", options.patch_weight, ...
                "sliver_correction", options.sliver_correction, ...
                "source_neumann", {options.source_neumann}, ...
                "source_dirichlet", {options.source_dirichlet}, ...
                "interface_law", options.interface_law, ...
                "slip_coefficient", options.slip_coefficient, ...
                "swap_transmission", options.swap_transmission, ...
                "source_slip", {options.source_slip});
            obj = InterfaceProblem([fluid_problem, solid_problem], interface, ...
                coupling, ...
                quadrature_patch = options.quadrature_patch, ...
                quadrature_path = options.quadrature_path, ...
                vertex_tol = options.vertex_tol, ...
                vertex_tol_relative = options.vertex_tol_relative, ...
                tip_mapping = options.tip_mapping, ...
                collapsed_normal = options.collapsed_normal, ...
                mapping_1 = options.mapping_1, mapping_2 = options.mapping_2);
        end

        function obj = stokes_biot(fluid_problem, biot_problem, interface, options)
            % Stokes-Biot transmission (quasi-static Biot, one time step):
            % side 1 = stress-velocity PEERS fluid (IncompressibleElasticity
            % Matrices layout, "pinned" closure: no ConstantBasis global
            % space), side 2 = Biot (BiotMatrices layout: local [sigma_p,
            % eta, gamma_p, q, p], one three-component trace
            % [eta_hat_x; eta_hat_y; p_hat] per face). The fluid pressure
            % level is fixed by R_p. Rows, with
            % nu_p the unit normal of Sigma out of the Biot side, nu_f =
            % -nu_p, t = (-nu_p2, nu_p1), cdt the effective time step and
            % h* = eta_hat* + Ext(eta*) the transported history:
            %   R_f (fluid trace rows, vector Legendre test): traction
            %       balance <Eh(sigma_f) + Eh(sigma_p), mu> = <T_1, mu>;
            %   R_n (eta_hat_x rows, scalar test): mass conservation
            %       <nu_p.(cdt Th(u_hat_f) - Th(eta_hat)) - cdt Eh(q_hat), l>
            %         = <-nu_p.h* + cdt nu_p.T_2, l>;
            %   R_t (eta_hat_y rows, scalar test): BJS friction on the fluid
            %       stress, <t.Eh(sigma_f) + (g/cdt) t.(cdt Th(u_hat_f)
            %       - Th(eta_hat)), l> = <-(g/cdt) t.h* + t.T_4, l>,
            %       g = friction_bjs;
            %   R_p (p_hat rows, scalar test): normal stress
            %       <Tp(p_hat) + nu_f.Eh(sigma_f), l> = <T_3, l>.
            % R_n and R_t are the rows of fsi's "navier" law with
            % slip_coefficient = friction_bjs/cdt, plus the Darcy flux in
            % R_n and the c_b p term (c_b = alpha/(2(lambda + mu))) in the
            % displacement path integral Th(eta_hat). Volume time terms
            % belong to the Biot problem's local matrices.
            %
            % Pore-pressure level: with s0 = 0 and alpha = 1 (to a relative
            % tolerance 1e-10 on s0 (lambda + mu) and |1 - alpha|) the shift
            % (p, p_hat, sigma_P, sigma_F) -> (p + c, p_hat + c,
            % sigma_P - c I, sigma_F - c I) leaves every equation
            % unchanged unless a boundary off Sigma fixes it: a prescribed
            % p_hat, or a natural displacement or velocity component c on a
            % face with n_c ~= 0 (slip and symmetry walls and rollers, with
            % n_c = 0, do not fix it). When no such face exists
            % (InterfaceProblem.solve decides it from the applied boundary
            % conditions), the coupling-owned scalar beta
            % (n_extra = 1) closes the level: the row
            %   |Omega_P,h|^-1 (p, 1)_{Omega_P,h} = pore_pressure_mean
            % and the column -cdt <beta, l> in R_n (beta a uniform normal
            % flux added to Eh(q_hat) on Sigma, the multiplier that
            % absorbs the discrete incompatibility of the data). Otherwise
            % beta is eliminated (beta = 0) and the system is the one
            % without it. solve reports beta as global_system.level_scalar.
            arguments
                fluid_problem (1,1) FEProblem
                biot_problem  (1,1) FEProblem
                interface AbstractCurve
                options.viscosity           (1,1) double = 1.0
                options.young_modulus       (1,1) double = 1.0
                options.poisson_ratio       (1,1) double = 0.25
                options.biot_alpha          (1,1) double = 1.0
                options.diffusion           (2,2) double = eye(2)
                % tau_p / (||kappa|| tau); empty = the biot_problem's
                % Data.stabilization_darcy (1 when absent)
                options.stabilization_darcy double {mustBeScalarOrEmpty} = []
                % tau_s / tau on the Biot side; empty = the biot_problem's
                % Data.stabilization_solid (1 when absent)
                options.stabilization_solid double {mustBeScalarOrEmpty} = []
                % s0, read only for the pore-pressure level; empty = the
                % biot_problem's Data.storage (0 when absent)
                options.storage double {mustBeScalarOrEmpty} = []
                % the mean of p over the Biot mesh, imposed only when the
                % pore-pressure level is free (see above)
                options.pore_pressure_mean  (1,1) double = 0.0
                options.friction_bjs        (1,1) double {mustBeNonnegative} = 0.0
                options.effective_time_step (1,1) double {mustBePositive, mustBeFinite} = 1.0   % cdt
                % the combinations w* of biot_step_history (fields
                % displacement_skeleton, stress, rotation, pressure; [], 0
                % or missing = zero), plus the optional field extension
                % (Ext(eta*) given directly, (2*d1) x n_patches, used
                % instead of the transport built from stress, pressure and
                % rotation); struct.empty = zero history, no transport
                options.history struct = struct.empty
                options.source_neumann   cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}       % T_1
                options.source_dirichlet cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}       % T_2
                options.source_pressure  (1,1) function_handle = @(x,y) 0*x.*y       % T_3
                options.source_slip      cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}       % T_4
                options.quadrature_patch double = []
                options.quadrature_path  double = []
                options.vertex_tol       (1,1) double = 1e-8
                options.vertex_tol_relative double = []
            end
            has_mu_c = any(arrayfun(@(sp) isa(sp.Basis, 'ConstantBasis'), ...
                fluid_problem.GlobalSpaces));
            assert(~has_mu_c, "InterfaceProblem:stokesBiotClosure", ...
                "stokes_biot needs the fluid's ""pinned"" closure (no ConstantBasis " + ...
                "global space): the fluid pressure level is fixed by the R_p row.");
            history_extension = [];
            if isfield(options.history, "extension")
                history_extension = options.history.extension;
            end
            stabilization_darcy = options.stabilization_darcy;
            if isempty(stabilization_darcy)
                stabilization_darcy = 1.0;
                if isstruct(biot_problem.Data) && isfield(biot_problem.Data, "stabilization_darcy")
                    stabilization_darcy = biot_problem.Data.stabilization_darcy;
                end
            end
            stabilization_solid = options.stabilization_solid;
            if isempty(stabilization_solid)
                stabilization_solid = 1.0;
                if isstruct(biot_problem.Data) && isfield(biot_problem.Data, "stabilization_solid")
                    stabilization_solid = biot_problem.Data.stabilization_solid;
                end
            end
            storage = options.storage;
            if isempty(storage)
                storage = 0;
                if isstruct(biot_problem.Data) && isfield(biot_problem.Data, "storage")
                    storage = biot_problem.Data.storage;
                end
            end
            coupling = struct("kind", "stokes_biot", "n_extra", 1, ...
                "viscosity", options.viscosity, ...
                "young_modulus", options.young_modulus, ...
                "poisson_ratio", options.poisson_ratio, ...
                "biot_alpha", options.biot_alpha, ...
                "diffusion", options.diffusion, ...
                "stabilization_darcy", stabilization_darcy, ...
                "stabilization_solid", stabilization_solid, ...
                "storage", storage, ...
                "pore_pressure_mean", options.pore_pressure_mean, ...
                "friction_bjs", options.friction_bjs, ...
                "effective_time_step", options.effective_time_step, ...
                "history_skeleton", history_field(options.history, "displacement_skeleton"), ...
                "history_stress", history_field(options.history, "stress"), ...
                "history_rotation", history_field(options.history, "rotation"), ...
                "history_pressure", history_field(options.history, "pressure"), ...
                "history_extension", history_extension, ...
                "source_neumann", {options.source_neumann}, ...
                "source_dirichlet", {options.source_dirichlet}, ...
                "source_pressure", options.source_pressure, ...
                "source_slip", {options.source_slip});
            obj = InterfaceProblem([fluid_problem, biot_problem], interface, ...
                coupling, ...
                quadrature_patch = options.quadrature_patch, ...
                quadrature_path = options.quadrature_path, ...
                vertex_tol = options.vertex_tol, ...
                vertex_tol_relative = options.vertex_tol_relative);
        end
    end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function h = history_field(history, name)
% a field of the fsi history option; [] or missing = 0
h = 0.0;
if isfield(history, name) && ~isempty(history.(name))
    h = history.(name);
end
end
