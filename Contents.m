% TPM-HDG  HDG methods with the Transfer Path Method for unfitted and dissimilar meshes
% Version 1.0.0 (R2024b)
%
% Setup
%   initpath                            - Add the toolbox folders to the path (run first).
%
% Problems and solutions
%   FEProblem                           - Variational problem: spaces, operator, boundary conditions.
%   FESolution                          - Solution of an FEProblem (fields and hybrid unknowns).
%   FESpace                             - Discrete space on a mesh, condensed or retained.
%   AffineOperator                      - Term list of integral forms and pre-assembled blocks.
%   BoundaryConditions                  - Weak, strong and transferred boundary conditions.
%   TransferEssentialBoundaryConditions - Dirichlet data transferred along TPM paths.
%   FEField                             - Discrete field; compute_l2_error, gradient, plot.
%   Datum                               - Data function in the form language.
%   Measure                             - Integration measure (dx, ds) of the form language.
%
% PDEs
%   DarcyProblem                        - Darcy flow.
%   StokesProblem                       - Stokes flow.
%   OseenProblem                        - Oseen flow (Stokes with convection).
%   NavierStokesProblem                 - Navier-Stokes flow (Picard iteration).
%   IncompressibleElasticityProblem     - Incompressible elasticity / Stokes, stress-velocity form (PEERS).
%   LinearElasticityProblem             - Linear elasticity, steady and time-dependent.
%   HyperelasticProblem                 - Hyperelastic solid, total Lagrangian, on the PEERS unknowns.
%   HyperelasticGradientProblem         - Hyperelastic solid, displacement gradient as mixed variable.
%   BiotProblem                         - Biot poroelasticity (PEERS total stress, Darcy flux).
%   hdg_solve_cartesian                 - Q^k HDG advection-diffusion on a CartesianGrid.
%   hdg_poisson_cartesian               - Q^k HDG Poisson on a CartesianGrid.
%
% Coupled problems (TPM)
%   InterfaceProblem                    - Two problems coupled across an interface.
%   FluidStructureCoupling              - Fluid-structure interaction (package).
%   StokesDarcyCoupling                 - Stokes-Darcy coupling, BJS friction (package).
%   StokesBiotCoupling                  - Stokes-Biot and Navier-Stokes-Biot coupling (package).
%
% Interface curves
%   AbstractCurve                       - Base class: partitions, transfer paths, mappings.
%   ParametrizedCurve                   - Curve given by a parametrization.
%   LevelSetCurve                       - Curve given implicitly by a level set.
%   CircularCurve                       - Circle.
%   AirfoilCurve                        - Joukowski airfoil.
%   KidneyCurve                         - Kidney-shaped level-set curve.
%   DeformedCurve                       - Image of a reference curve under a displacement field.
%   plot_transfer_paths                 - Plot the transfer paths of a partition.
%
% Meshes
%   AbstractGrid                        - Base class of the grids.
%   SimplicialGrid                      - Triangular mesh with boundary tags.
%   StructuredSimplicialGrid            - Structured triangulation of a rectangle.
%   CartesianGrid                       - Tensor-product grid of rectangles.
%   interface_meshes                    - Meshes on both sides of an interface curve.
%   interface_polygon                   - A mesh's own interface facets, as a PolygonalCurve.
%   immerse_mesh                        - Carve the domain {phi < 0} out of a background grid.
%
% Bases
%   PEERSBasis                          - PEERS matrix-valued stress basis (Dubiner + bubbles).
%   QkBasis                             - Tensor-product Legendre basis of Q^k.
%   LagrangeBasis                       - Continuous Lagrange basis (linear only).
%
% Verification
%   run_suite                           - Run every test in tests/.
