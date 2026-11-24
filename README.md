# TPM-HDG: Hybridizable Discontinuous Galerkin methods for unfitted and dissimilar meshes

A MATLAB toolbox for 2D Hybridizable Discontinuous Galerkin (HDG) methods,
with the Transfer Path Method (TPM) of Cockburn and Solano for unfitted
boundaries and for coupling problems posed on dissimilar, non-matching
meshes. It covers Darcy, Stokes, Oseen/Navier–Stokes, linear elasticity
(compressible and incompressible), hyperelastic solids and Biot
poroelasticity, steady and time-dependent; the Stokes–Darcy,
fluid–structure, Stokes–Biot and Navier–Stokes–Biot couplings; and HDG
with tensor-product (Q^k) bases on Cartesian grids.

The HDG machinery follows Fu, Gatica and Sayas (2015); the implementation
of TPM started in 2020 during my MSc thesis.

## Requirements and setup

- MATLAB R2024b or later, with the Symbolic Math and Parallel Computing toolboxes.
- The mesh generator [`mesh2d`](https://github.com/jaimemanriquezr/mesh2d)
  (a modified copy of D. Engwirda's MESH2D), a git submodule. Clone with
  `git clone --recurse-submodules https://github.com/jaimemanriquezr/TPM-HDG.git`,
  or run `git submodule update --init` in an existing clone. The `.mltbx`
  file already contains it.

TPM-HDG is a MATLAB package (`resources/mpackage.json`, name `tpmhdg`,
display name TPM-HDG). It installs in one of three ways; each puts `src/`,
`examples/` and `tests/` (with their subfolders) on the path.

1. **With the MATLAB package manager.** From a clone (with the submodule):
   `mpminstall("path/to/TPM-HDG")`. Uninstall with `mpmuninstall("tpmhdg")`.
2. **From source.** In a clone, run `initpath` once per session from the
   repository root; nothing is installed.
3. **The toolbox file from Releases.** Download `TPM-HDG.mltbx` from the
   repository's Releases page and open it in MATLAB, or run
   `matlab.addons.install("TPM-HDG.mltbx")`. It includes `mesh2d`.
   Uninstall from the Add-On Manager or with
   `matlab.addons.uninstall("TPM-HDG")`.

`Contents.m` lists the main classes and functions by area.

### Development

Open the MATLAB project `TPM-HDG.prj` (`openProject("path/to/TPM-HDG")`);
the project sets the path, and closing it removes it. MATLAB warns on
opening that the package root is on the path but the package is not
installed; the warning is harmless. Do not also install the package
(route 1) while working in the project.

`buildtool package` (see `buildfile.m`) checks the package definition
against the tree and writes `release/TPM-HDG.mltbx`, the file attached to
each release; `buildtool folders` updates the package's folder list after
folders are added or removed.

## Layout

| Directory | Contents |
| --- | --- |
| `src/fem/` | Finite element core: bases (including the tensor-product `QkBasis`), quadrature, the form language, `FEProblem`/`FESolution`, matrix builders, boundary conditions, the PDE classes (`DarcyProblem`, `StokesProblem`, `OseenProblem`, `LinearElasticityProblem`, `IncompressibleElasticityProblem`, `NavierStokesProblem`, `BiotProblem`, and the hyperelastic solids `HyperelasticProblem` (PEERS, Newton, Kelvin–Voigt viscosity) and `HyperelasticGradientProblem`), `InterfaceProblem`, and the Q^k HDG solvers on Cartesian grids (`hdg_solve_cartesian`, `hdg_poisson_cartesian`) |
| `src/dissimilar-coupling/` | TPM: interface curves (parametrized and level-set; `DeformedCurve`, the current interface of a Lagrangian solid), induced partitions of two meshes, transmission conditions, and the `StokesDarcyCoupling`, `FluidStructureCoupling` and `StokesBiotCoupling` solvers; `InterfaceProblem.fsi` also takes a per-side curve pair {current, reference} sharing the parameter, the fluid traction scaled by the Nanson factor |
| `src/mesh-generation/` | `SimplicialGrid` and `CartesianGrid` (both `AbstractGrid`s) and mesh constructors: structured grids, fitted splits, immersed (carved) meshes, `interface_meshes`; `mesh2d` (submodule) |
| `examples/` | Documented, runnable drivers (below) |
| `tests/` | Verification scripts; `run_suite` runs them all |
| `resources/` | The package definition `mpackage.json` and the MATLAB project files |
| `buildfile.m`, `Contents.m`, `initpath.m` | Build tasks, the function index, and the path setup for use from source |

## Two ways to write a problem

A problem is an `FEProblem`: a set of `FESpace`s (each declaring whether it
is statically condensed), an operator, and boundary conditions. `solve`
dispatches on content.

1. **Declarative.** Write the bilinear and linear forms in the form
   language (`get_trial_functions`, `Measure`, `dot`, `divergence`, ...);
   the toolbox assembles, condenses, applies the boundary conditions and
   recovers the local fields.
2. **PDE classes.** `DarcyProblem`, `StokesProblem`, etc. provide the
   matrices directly, including their time discretizations
   (Newmark, one-step, theta schemes).

Both paths return an `FESolution` whose `Unknowns.field` / `Unknowns.hybrid`
hold `FEField`s; errors are computed with `compute_l2_error` against
exact solutions.

## Examples

### A PDE class: Stokes on the unit square

```matlab
k = 2;
data = StokesProblem.get_manufactured_solution("x^2*y^2", "x*y", 1.0);
problem = StokesProblem(viscosity = data.nu, source_term = data.f, datum_dirichlet = data.u);
mesh = StructuredSimplicialGrid(20, 20);
solution = solve(problem, mesh, Order = k);
[e, e_hat] = solution.compute_l2_error({data.L, data.u, data.p}, {data.u});
e_u_hat = e_hat(1);  % e_hat(2) is the pressure-mean multiplier, compared against zero
```

The second hybrid unknown of `StokesProblem` is the element-wise
pressure-mean multiplier, which has no exact datum here, so `e_hat(2)` is
its norm, not an error.

### The declarative path: Darcy in the form language

```matlab
k = 2;
g = Datum(@(x, y) sin(pi*x).*sin(pi*y));
f = Datum(@(x, y) 2*pi^2*sin(pi*x).*sin(pi*y));
Pk = DubinerBasis(k); VPk = Pk.to_vector(); Mk = LegendreBasis(k);
mesh = StructuredSimplicialGrid(16, 16).add_stabilization(1);

problem = FEProblem(Spaces = [FESpace([VPk, Pk], Condense = true), FESpace(Mk)]);
problem.BoundaryConditions = WeakEssentialBoundaryConditions(g, "Dirichlet");
[q, u, u_hat] = get_trial_functions(problem.Space);
[r, w, v_hat] = get_test_functions(problem.Space);
dx = Measure(preset = "dx");  dS = Measure(preset = "dS");  n = ExteriorNormal();

problem.Operator = AffineOperator( ...
      dot(q, r)*dx - u*divergence(r)*dx + u_hat*dot(r, n)*dS ...
    + divergence(q)*w*dx + u*w*dS - u_hat*w*dS + f*w*dx ...
    - (dot(q, n)*v_hat*dS + u*v_hat*dS - u_hat*v_hat*dS));
solution = solve(problem, mesh, k);
e_u = compute_l2_error(solution.Unknowns.field(2), mesh, g);
```

### An unfitted boundary: Darcy on a level-set domain

The domain `{phi < 0}` of a `KidneyCurve` is carved out of a background
grid; the Dirichlet datum is transferred from the true boundary to the
polygonal computational boundary along transfer paths. See
`examples/hdg/pde_darcy_levelset_example.m`.

```matlab
kidney = KidneyCurve();
mesh = immerse_mesh(StructuredSimplicialGrid(32, 32, -0.8, 0.8, -0.8, 0.8), kidney);
problem.BoundaryConditions = TransferEssentialBoundaryConditions.darcy(u_exact, "Dirichlet", kidney, kappa);
```

### Coupling on dissimilar meshes: fluid–structure interaction

Stokes flow coupled with linear elasticity across a flat interface, the
two meshes generated independently:

```matlab
k = 1;  h = 1/4;
interface = FlatCurveHorizontal(0.5);
[fluid_mesh, solid_mesh] = interface_meshes(interface, "fitted", h);
[fluid_problem, solid_problem, sources, exact] = FluidStructureCoupling.load_problems( ...
    k, interface, 0, 1e-3, "x*y", "x^3 + y^3", ["x"; "y"], 0, 0.4, 0, 1.0, 0.25, 0);
[fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
    fluid_problem, solid_problem, fluid_mesh, solid_mesh, interface, sources{:});
e_fluid = fluid_solution.compute_l2_error({exact.sigma_F, exact.u_F, exact.gamma_F}, {exact.u_F});
```

`interface_meshes` accepts the families `"fitted"`, `"interior"`,
`"protrude_fluid"` and `"protrude_solid"`, which place the two meshes
against the interface with or without a gap. `FluidStructureCoupling.solve`
takes a scheme (`"steady"`, `"onestep"`, `"newmark"`) and an
`interface_law` (`"noslip"`, `"navier"`, `"full"`).

With scheme `"steady"`, `time_step = Inf` is the stationary problem
(kinematic row <u_F, mu> = <T2, mu>, no solid term) and requires
`pressure_mean = "constrained"`; see `tests/test_interface_fsi_stationary.m`.

### Index of examples

| File | Problem |
|---|---|
| `hdg/pde_{darcy,stokes}_example.m` | One solve of each PDE class on the unit square |
| `hdg/pde_elasticity_example.m` | Time-dependent linear elasticity on the unit square |
| `hdg/pde_{darcy,stokes,elasticity}_convergence.m` | Convergence tables (EOC ≈ k+1) |
| `hdg/pde_elasticity_mixed_bc.m` | Elasticity with mixed Dirichlet/Neumann boundaries |
| `hdg/pde_stokes_2_example.m` | Time-dependent Stokes |
| `hdg/pde_stokes_2_convergence.m` | Steady Stokes, stress–velocity–vorticity form |
| `hdg/dsl_darcy_example.m` | Darcy written in the form language |
| `hdg/dsl_stokes_2_convergence.m` | Stress–velocity–vorticity Stokes in the form language (declarative path) |
| `hdg/dsl_darcy_qk_example.m` | Darcy with Q^k on a `CartesianGrid`, in the form language |
| `hdg/qk_vs_pk_poisson.m` | Poisson with Q^k on a `CartesianGrid` against P^k on triangles |
| `hdg/pde_darcy_levelset_example.m` | Unfitted Darcy on a level-set domain |
| `hdg/pde_navier_stokes_convergence.m` | Navier–Stokes by Picard iteration |
| `hdg/pde_navier_stokes_unfitted_convergence.m` | Navier–Stokes on an immersed disk |
| `hdg/cylinder_benchmark.m` | Schäfer–Turek flow past a cylinder, unfitted |
| `meshes/*.m` | Mesh couplings: fitted, unfitted, closed interfaces |
| `tpm/tpm_meshes_example.m` | Basic domains: unit square, L-shape, square with a hole |
| `tpm/stokes_darcy_flat_{example,convergence}.m` | Stokes–Darcy with Beavers–Joseph–Saffman |
| `tpm/FSI_{flat,smooth,closed}_{example,convergence}.m` | FSI across flat, curved and closed interfaces |
| `tpm/FSI_closed_unfitted_convergence.m` | Convergence of FSI with an unfitted closed interface |
| `tpm/FSI_interface_example.m` | The FSI solve written out on `InterfaceProblem` |
| `tpm/FSI_navier_slip_example.m` | FSI with a Navier slip law |
| `tpm/FSI_dsl_example.m` | FSI with both subproblems in the form language on `InterfaceProblem` |
| `tpm/FSI_dsl_flat_convergence.m`, `tpm/FSI_dsl_native_flat_convergence.m` | FSI convergence on a flat interface with a gap closing as h^s (form language with `dot` / with MATLAB-native products) |
| `tpm/FSI_turek_hron_benchmark.m` | Turek–Hron FSI flag benchmark |
| `tpm/dsl-fsi/FSI_dsl_{curved,gap,lambda}_convergence.m` | FSI convergence on a curved interface / with a gap h^s / towards the incompressible solid limit |
| `tpm/dsl-fsi/fsi_dsl_picard_convergence.m` | Navier–Stokes/elasticity FSI by Picard iteration |

`tpm/dsl-fsi/fsi_dsl_*.m` are helper functions of these drivers, not drivers.

### Plotting

`plot(mesh)` / `plot(mesh, ax)` draws a `SimplicialGrid`;
`plot(curve, ax, Color = ...)` any curve; `plot(field, mesh)` an `FEField`;
`plot(solution)` a whole `FESolution` (all unknowns in a tiled layout,
`components = "all" | "tiled" | "first"`); `plot(datum, mesh)` data; and
`plot_transfer_paths(ax, interface, mesh, partition)` the TPM transfer
paths of an induced partition (`get_partitions`). `help <class>/plot`
lists the options.

## Tests

`run_suite` runs every `tests/test_*.m` and reports those that throw or
print `N check(s) failed` with N > 0; each script also runs on its own.
It assumes the toolbox is on the path (installed, the project open, or
`initpath` run); a healthy run reports `SUITE: 0/N failed`. Patch tests
require the exact solution to be reproduced to round-off: relative error
at most 1e-10, or 100 times the round-off of the global solve when that is
larger (measured by one step of iterative refinement,
`tests/helpers/solve_roundoff.m`; with a gap of order h at k = 3 the
transfer systems reach condition numbers of 1e9 to 1e10). In `test_transfer_oseen_curved` the
stress-form rates of sigma and gamma are reported rather than asserted: its
carved meshes leave transfer paths of order h, outside the regime in which
the method's rates are guaranteed.

## References

1. Z. Fu, L. F. Gatica, F.-J. Sayas. [Algorithm 949: MATLAB Tools for HDG in Three Dimensions](https://doi.org/10.1145/2658992). ACM Trans. Math. Softw., 2015.
2. M. Solano, S. Terrana, N.-C. Nguyen, J. Peraire. [An HDG method for dissimilar meshes](https://doi.org/10.1093/imanum/drab059). IMA J. Numer. Anal., 2022.
3. JM, N.-C. Nguyen, M. Solano. [A dissimilar non-matching HDG discretization for Stokes flows](https://doi.org/10.1016/j.cma.2022.115292)\*. Comput. Methods Appl. Mech. Engrg., 2022.
4. I. Bermúdez, JM, M. Solano. [A hybridizable discontinuous Galerkin method for Stokes/Darcy coupling on dissimilar meshes](https://doi.org/10.1093/imanum/drae109)\*. IMA J. Numer. Anal., 2025.

\*Results obtained with this code.
