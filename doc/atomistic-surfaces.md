# Fixed-site atomistic surfaces

Apothesis can now construct a slab as explicit three-dimensional crystal sites.
Each site has a stable ID, Cartesian position, species, occupancy, fixed-substrate
flag and geometrical neighbors. Vacancies remain candidate sites, so an interior
vacancy does not collapse a column of atoms above it.

The fixed-site atomistic backend supports adsorption, diffusion and desorption
KMC events. Coordinates and the candidate-site graph remain fixed; events change
individual occupancies and species. The height-field backend remains separate.
There is no molecular dynamics, geometry relaxation or atomistic reaction backend.

## Build a copper slab

```sh
./apothesis --validate-input examples/atomistic-cu.json
./apothesis --build-surface examples/atomistic-cu.json --output copper.xyz
```

The example produces 576 Cu atoms on an FCC slab containing 6 × 6 × 4 conventional
unit cells. Lengths are in angstrom. The bottom unit cell in z (144 atoms) is marked
fixed. A total vacuum thickness of 15 angstrom is divided equally below and above
the repeated unit-cell volume. The cell is periodic in x and y, open in z.
Without `--output`, the filename is `surface.xyz` in the working directory.

Extended XYZ output contains the cell vectors, periodic flags, element, position,
stable site ID and fixed flag. Vacant sites are omitted from the atom records;
their IDs and topology remain present in memory. The file is a geometry export,
not a restart file for the complete candidate-site graph.

## Input

```json
{
  "schema_version": 1,
  "lattice": {
    "representation": "atomistic",
    "type": "FCC",
    "parameters": {
      "lattice_constant": 3.615,
      "repetitions": [6, 6, 4],
      "initial_species": "Cu",
      "periodic": [true, true, false],
      "vacuum": 15,
      "fixed_unit_cells": 1,
      "vacancies": [100]
    }
  },
  "species": {"Cu": {"parameters": {"element": "Cu"}}}
}
```

Geometry-only input needs `schema_version`, `lattice` and `species`.
For KMC, supply all four additional sections: `system`, `simulation`, `processes`
and `output`. If any is supplied, the parser requires and validates all four.
`--build-surface` still only exports the initial geometry.
Additional parameter values remain in the parsed configuration document.
The initial geometry uses a single atomic species. An explicit chemical `element`
is required for that species; a molecular formula is not interpreted as an atom.

| Parameter | Meaning |
| --- | --- |
| `lattice_constant` | Positive lattice constant a, in angstrom. |
| `repetitions` | Three positive integers: numbers of unit cells along each lattice vector. |
| `initial_species` | A declared species with a chemical element symbol. |
| `length_unit` | Optional `angstrom`; other units are rejected. |
| `orientation` | Optional `100` for cubic structures, `0001` for HCP. |
| `periodic` | Three booleans; default `[true, true, false]`. |
| `origin` | Optional Cartesian translation `[x, y, z]`, in angstrom; negative values are allowed. |
| `vacuum` | Nonnegative total extra cell length along z; default zero. Requires nonperiodic z. |
| `occupied_unit_cells` | Initially fill the bottom nz unit cells; default fills all generated cells. Cells above this value are empty growth sites. |
| `fixed_unit_cells` | Number of bottom unit cells marked fixed, from zero through nz. |
| `vacancies` | Unique, zero-based site IDs to leave empty. |
| `neighbor_cutoff` | Optional positive geometrical cutoff in angstrom. |
| `c_over_a` | HCP c/a ratio, default sqrt(8/3). |

Supported bases: SimpleCubic (1 site/cell), BCC (2), FCC (4), Diamond (8), and
HCP (2). Cubic surfaces expose the {100} family with surface normal along the
Cartesian z axis. HCP uses the basal plane with nonorthogonal in-plane vectors.
Arbitrary Miller cuts, XYZ import, mixed-species bases and amorphous structures
are not part of this first geometry implementation. Height-field parameters such
as `initial_height` are rejected for atomistic input.

Site IDs are assigned with basis index varying fastest, then x, y and z unit-cell
indices. This ordering is independent of which sites are occupied. Removing a
site via `vacancies` changes atom count and occupied coordination, not site count.

## Neighbors and boundaries

Neighbors are computed geometrically with a spatial bin index. Periodic distances
use adjacent-image checks, including the skew HCP cell. Neighbor links are unique
site IDs and symmetric; periodic cell heights must exceed twice the cutoff to
avoid multiple-image ambiguity. Cells are limited to one million candidate sites.

The default cutoff is 1.01 times the ideal first-neighbor distance. For HCP it
includes both the basal and interplane distances; a nonideal c/a ratio may require
an explicit cutoff to select the intended shell. Geometric coordination includes
vacant candidate sites; `occupiedCoordination(id)` counts occupied neighbors.
There is no automatic claim that every undercoordinated atom is an adsorption site.

## Run atomistic KMC

```sh
./apothesis --validate-input examples/atomistic-kmc.json
./apothesis --input examples/atomistic-kmc.json
```

The example has 256 FCC candidate sites, 128 initial Cu atoms, 64 fixed substrate
atoms and 128 initially empty growth sites. It includes all three processes.
The rates are illustrative; they are not calibrated copper kinetics.

`occupied_unit_cells` reserves empty crystal cells above the substrate. `vacuum`
adds empty simulation-box space only, not candidate sites. Growth is bounded by
the generated graph; add enough repetitions in z for the intended deposition.
A fully occupied slab has no initial adsorption events until a vacancy is exposed.

### Event eligibility

The current gas reservoir is on the **+z side**, so KMC requires nonperiodic z.
Empty candidate sites in the topmost plane seed a flood fill through neighboring
empty sites. This defines gas-connected vacancies. An occupied atom is exposed
if it is in the topmost candidate plane or neighbors a gas-connected vacancy.
This is a vacancy-graph accessibility model, not a ballistic collision model.

- **Adsorption**: fills one nonfixed, gas-connected empty site with the product
  species. By default it requires at least one occupied supporting neighbor.
- **Desorption**: empties one exposed, nonfixed site matching the reactant species.
  Buried atoms and atoms beside sealed internal vacancies are not exposed.
- **Diffusion**: moves one exposed, nonfixed atom to a neighboring nonfixed,
  gas-connected vacancy. The destination must retain at least one supporting
  occupied neighbor after the source is removed. Atom count and species are
  conserved, and the source becomes vacant.

Every occupied atom must have at least **one occupied first-nearest neighbor**.
An initial state containing isolated atoms is rejected. Adsorption and diffusion
cannot set `min_coordination` below one. Desorption and diffusion are rejected if
removing the source would leave any remaining atom with zero bonds. A hop is
allowed when its destination simultaneously restores the dependent atom's bond.
This rule applies to all species and fixed substrate atoms as well as the film.
When occupied fixed substrate sites are present, every occupied atom must also
be connected to at least one of them through occupied first-neighbor bonds.
Initial floating clusters are rejected, and events that detach an entire bonded
cluster are blocked. A simultaneous hop may replace the connection. If no fixed
substrate is specified, only the local one-bond minimum is enforced.

Fixed sites cannot change occupancy. The geometric neighbor graph remains fixed;
occupied coordination and the event catalog are rebuilt after each event. This
initial implementation prioritizes correctness over large-system throughput.

Process definitions use the same structured reactants/products as other inputs.
Only single-atom processes with coefficients of one are currently supported.
Adsorption/desorption gas and surface species must declare the same element;
different labels for the same element are allowed. The initial substrate may
have a different element from the deposited atoms, and XYZ exports each element.
Molecular dissociation, chemical reactions and multi-site moves are rejected.

Custom process `parameters` may specify `min_coordination` and `max_coordination`.
For adsorption these apply to the target; for diffusion to the target after the
source is removed; for desorption to the source. Defaults are minimum one for
adsorption/diffusion, zero for desorption, and no practical maximum. Setting the
minimum below one is rejected for adsorption and diffusion. The previous
`all_neighbors` flag must be false or absent: each eligible site or directed hop
is already represented as a separate event.

### Rates, time and stopping

Supported rate models:

- `constant`: nonnegative `value` in **s^-1 per eligible site**, or **per directed
  source–destination pair** for diffusion. There are no hidden neighbor factors
  or conversions from ML/s.
- `arrhenius`: nonnegative `prefactor` in s^-1 and `activation_energy` in J/mol;
  optional nonnegative `bond_energy` in J/mol per occupied neighbor. The rate is
  `prefactor * exp(-(activation_energy + bond_energy*n)/(R*T))`, with T in K.
  For desorption/diffusion, n is the source's pre-event occupied coordination;
  for adsorption, n is the target's coordination. This is an explicit optional
  bond-counting model, not a calibrated energy landscape. For a hop, put the
  migration barrier in `activation_energy`; `migration_energy` is rejected.

The `simple` gas-flux model is not supported by this backend yet. Model parameters
must be supplied and validated for the physical system being simulated.
`system.parameters.temperature` is required; pressure is optional and does not
alter constant or Arrhenius rates. Extra custom parameter objects are preserved.

The direct KMC algorithm selects individual events with probability rate/total
rate and samples `dt = -log(u)/total_rate` using the **pre-event** total rate.
Events beyond `start_time + duration` are not executed. With no eligible nonzero
rates, the state stays unchanged through the remaining time and stops cleanly.
See the [surface KMC practical guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC6465329/)
for the underlying direct-method formulation.

`simulation.random_seed` supplies reproducible runs within this implementation;
zero/absent chooses an automatically generated seed, recorded in the log.
Optional `simulation.max_events` stops after that many events; zero/absent means
no event-count limit. A capped run reports `event_limit` and its actual elapsed
time instead of pretending to have reached the time horizon.

### Outputs

Files are written in the working directory:

- `atomistic.log`: seed, actual simulation times, atom count, total event count,
  per-process counters and stop reason (`end_time`, `no_events` or `event_limit`).
- `atomistic_initial.xyz` and `atomistic_final.xyz`: initial and final atoms with
  cell, periodic flags, site IDs, fixed flags and time metadata.
- `atomistic_trajectory.xyz`: concatenated extended XYZ frames for visualization.

`output.log_interval` and `output.lattice_interval` are positive intervals in
seconds. Intermediate output is sampled at the first event after an interval has
elapsed, labeled with the actual event time; the initial and final states are
always included. Height-field coverage reporting is not supported here.
The exported trajectory is not a complete restart of vacancies, topology and RNG.
