# CrystalSurface: first implementation

Base: nixeimar/Apothesis commit 373e53e7523a55c4fe472e71f98f810e7b052918.
License: GPL-3.0-or-later, as for Apothesis.

This C++17 module represents adsorption candidates on ideal, unreconstructed
FCC(111), FCC(100), and HCP(0001) slabs. It is built as apothesis_crystal by
CMake and linked to Apothesis. qmake includes the source as well.

## Scope and integration status

The module and a separate working direct-method KMC example are implemented.
The existing input.kmc parser and legacy height-column process dispatch do NOT
use CrystalSurface yet. Linking the library does not change existing runs.
Do not select this geometry through the old FCC/HCP input options expecting
top/bridge/hollow kinetics. This is the first fixed-slab module, not a completed
replacement of the deposition engine.

Supported:
- explicit slab atoms, periodic cell vectors, top/bridge/hollow candidates;
- image-aware supporting atoms and directed transition edges;
- user-enabled metastable site types (all disabled initially);
- configurable hard-core exclusions, adsorption, desorption and diffusion;
- conservative complete event invalidation after changes;
- atom/site/transition CSV export;
- deterministic example with adsorption, desorption and hopping.

Not supported: arbitrary Miller indices, reconstructed or mixed-material
surfaces, relaxation, species-specific footprints, automatic energy barriers,
multilayer growth, or updating the substrate during deposition.

## Geometry conventions

Coordinates are a local surface frame, with outward normal +z and top layer z=0.
FCC a is the conventional cubic lattice parameter. HCP a and c are independent.
The square FCC(100) surface spacing is a/sqrt(2); its layers alternate offsets
(0,0) and (1/2,1/2), with vertical separation a/2.
FCC(111) primitive vectors enclose 60 degrees, spacing a/sqrt(2), layer
separation a/sqrt(3), and downward stacking offsets 0, 1/3, 2/3 along both
primitive vectors. HCP(0001) has spacing a, separation c/2 and offsets 0,1/3.
At least three slab layers are required.

Per primitive cell FCC(100) has 1 top, 2 bridge, 1 fourfold candidate.
The triangular surfaces have 1 top, 3 bridge, 1 hcp and 1 fcc hollow.
The hcp hollow lies above the second substrate layer; for FCC(111) the fcc
hollow lies above the third. On HCP the fcc label denotes the alternative
stacking continuation, not a substrate atom in the third layer.
Positions are projected geometric centres at z=0, not relaxed adsorption
heights. In a 1x1 cell multiple support images can refer to the same atom ID:
the image shifts must be retained.

## Physical model is explicit

enable(type) declares a candidate to be an allowed occupied state. Bridges
remain disabled in the example; geometry alone cannot decide whether a bridge
is a metastable state or a saddle point.

setExclusionDistance(r) imposes a user-selected species-independent lateral
hard-core distance. This is an approximation, not an inferred physical law.
The threshold is strict (<r), must be nonnegative and below half the shortest
cell vector, and cannot be changed with adsorbates present. Initially no
cross-site exclusion is imposed.

connect(type1,type2,cutoff) registers geometric transition candidates in both
directions, with cutoff <= nearest substrate-atom spacing. It does not infer
a barrier or prove an allowed pathway. Periodic image edges are distinct;
repeated calls do not duplicate them. Self-image edges in very small cells
are retained geometrically but cannot change occupancy and are not KMC hops.

For a calibrated process model, associate each directed transition with its
species and local-environment-dependent prefactor/barrier. Use consistent
energy units and ensure reverse rates satisfy detailed balance where
appropriate. Never use the example's dimensionless rates as material data.

Change.affected contains ALL site IDs intentionally: recomputing every event
is safe even for rates with nonlocal dependence. The reference implementation
uses quadratic geometric searches and is intended for modest cells. A scalable
version should use spatial bins and explicit rate-dependency neighbourhoods.

## Building and running

From the repository root:

    cmake -S . -B build-crystal -DCMAKE_BUILD_TYPE=Debug
    cmake --build build-crystal -j
    ctest --test-dir build-crystal --output-on-failure
    ./build-crystal/crystal_surface_demo fcc111 fcc111

Other faces are fcc100 and hcp0001. The optional second argument is the output
prefix. The example writes PREFIX_atoms.csv, PREFIX_sites.csv and
PREFIX_edges.csv (initial geometry), then executes 1000 toy KMC events.
It uses seed 12345, adsorption 1, desorption 0.2, hopping 2, in dimensionless
units, and verifies exclusions at the end. These are smoke simulations,
not predictions for a particular material.

Without CMake:

    g++ -std=c++17 -Isrc/crystal src/crystal/crystal_surface.cpp tests/crystal_surface_tests.cpp -o crystal_tests
    ./crystal_tests
    g++ -std=c++17 -Isrc/crystal src/crystal/crystal_surface.cpp examples/crystal_surface_demo.cpp -o crystal_demo
    ./crystal_demo fcc111

## Next integration step

Add an explicit crystal-mode input section and a graph-event scheduler that
uses CrystalSurface occupancy and transitions. Avoid adapting Site::height
to mean a hollow/bridge position: its old growth semantics are incompatible.
The example demonstrates the complete recompute-after-change scheduler.
Before production use, supply a calibrated material-specific event catalogue,
validate equilibrium/kinetics against a reference, and define film growth
rules that create or remove adsorption sites. The substrate is fixed here.
