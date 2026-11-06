# Geometry-based surface events

The `surface_kmc_3d` engine now has a geometry-based matcher, selected by
`model.parameters.surface_sites: "gas_exposed"`. This is the default for new
inputs that omit `surface_sites`. It uses actual fixed 3D coordinates and a
neighbor graph, without imposing a square surface or a vertical column per site.

The existing ZGB example explicitly selects `column_tops`. That mode preserves
the paper's solid-on-solid restrictions, same-height square-grid neighbors and
optimized incremental catalogs. Changing to `gas_exposed` changes the available
surface sites and hence the physical model; it is not an equivalent ZGB run.

## Geometry

Generated slabs support `SimpleCubic`, `BCC`, `FCC`, `HCP` and `Diamond`, using the
existing `lattice_constant`, `repetitions`, `occupied_unit_cells`, fixed layers,
vacancies, origin and lateral periodicity settings. z must remain open to the
+z gas reservoir. FCC/BCC and diamond use their conventional-cell bases; HCP
retains its skew in-plane cell and basal orientation. Other cuts can be supplied
as explicit candidate coordinates.

Run an [FCC example with a three-site reaction](../examples/surface-fcc-multisite.json):

```sh
./apothesis --input examples/surface-fcc-multisite.json
```

A custom geometry supplies an orthorhombic box and all candidate positions:

```json
"lattice": {
  "representation": "atomistic",
  "type": "Custom",
  "parameters": {
    "cell": [8, 4, 4],
    "periodic": [false, false, false],
    "neighbor_cutoff": 1.01,
    "site_volume": 1,
    "initial_species": "substrate",
    "sites": [
      {"position": [0, 1, 0], "occupied": true, "fixed": true},
      {"position": [1, 1, 0], "occupied": true, "fixed": true},
      {"position": [1, 1, 1], "occupied": false}
    ]
  }
}
```

Coordinates and box lengths are in angstrom; `site_volume` is in cubic angstrom
per effective solid unit. An individual site may override `species`. Site IDs
are the array indices and remain stable when occupancy changes. Vacant growth
positions must be supplied too. Custom positions stay inside the box; periodic
upper faces are excluded to avoid duplicate images. Arbitrary fixed coordinate
networks and mixed initial species are allowed. Custom cells are orthorhombic;
this interface does not import XYZ files or create off-lattice trajectories.

The complete [custom-coordinate example](../examples/surface-custom-coordinates.json)
contains a supported overhang and solid desorption. `--build-surface` also accepts
custom inputs. All models keep finite candidate positions; coordinates do not
move or relax during KMC.

## Exposure and connectivity

The **support graph** uses the geometry's first-neighbor cutoff. Empty candidate
sites at the top candidate plane seed a flood fill through empty neighbors.
Adsorbates occupy empty solid nodes and do not block this vacancy-graph gas path.
This defines access from +z; it is a topological access model, not a ballistic
collision or molecular diffusion calculation.

A vacant binding site must be gas-connected and have an occupied solid first
neighbor. A solid reactant must be nonfixed and gas-exposed: adjacent to that
empty gas region, or on the top candidate plane. Buried cavities can exist, but
sealed vacancies cannot adsorb. Supported overhangs and steps are allowed.

Each occupied solid must retain at least one occupied first neighbor. If any
fixed occupied substrate exists, every solid must remain connected to it.
Simultaneous updates are checked as a whole, so a diffusion move can replace a
supporting bridge, but cannot detach a dimer that still has internal bonds.
Existing and newly produced adsorbates must remain gas-accessible and supported.
Fixed substrate sites cannot be transformed or removed. Failed transitions are
excluded from the event propensity before the KMC clock is sampled.

## Connected multi-site patterns

Slots are local event roles, with contiguous indices starting at zero. All slots
must appear exactly once on each side of the process. A site cannot fill two
slots in the same event. Gas terms remain reservoir labels without site indices.
The `links` array states the required geometrical relations:

```json
"parameters": {
  "counting": "unique_site_sets",
  "links": [
    {"sites": [0, 1], "neighbors": "first", "same_height": false},
    {"sites": [1, 2], "neighbors": "first", "same_height": false}
  ]
}
```

This matches a connected three-site chain, with no required relation between
slots 0 and 2. Add a third link to require a triangle. Multiple distance links
can constrain angles indirectly through the corresponding pair distances. The
matcher checks all links, including cycle-closing links, on the actual geometry.
It supports up to 64 slots per pattern as a resource limit, rather than the old
two-slot restriction. Disconnected patterns are rejected.

Each link supports:

- `neighbors: "first"`: the geometry's support-neighbor graph, with arbitrary
  coordination number and no automatic height restriction.
- `neighbors: "lateral_first"`: that graph filtered to equal z.
- `neighbors: "distance"`: a separate geometric distance range in angstrom,
  requiring `max_distance`, with optional `min_distance` (default zero).
- `same_height: true`: additionally require equal Cartesian z. False/default
  permits unequal heights. `lateral_first` already implies equal height.

For example, simple-cubic second neighbors at lattice constant 1 can be selected
with `neighbors: "distance", min_distance: 1.4, max_distance: 1.43`. This does not
change the first-neighbor graph used to judge solid support or gas access.
Distances include periodic images and the generated HCP skew cell. Periodic
cell heights must exceed twice the largest requested distance; ambiguous
multiple-image interactions are rejected rather than counted as duplicate atoms.

For a two-slot rule, `neighbors` and `same_height` can be written directly in
process `parameters` as shorthand for one link between slots 0 and 1. Do not mix
shorthand and an explicit `links` array.

## Reactants and simultaneous updates

A `surface` reactant matches an adsorbate on a binding node; a `solid` reactant
matches an exposed occupied solid node; a vacancy matches an empty binding node.
A `solid` product occupies the mapped crystal node, a `surface` product creates
an adsorbate there, and a vacancy clears the mapped state. A process can create,
remove or transform multiple solid units simultaneously, subject to the support
and access checks. Newly accessible sites are recomputed from the resulting
occupancy instead of exposing only the node directly above a deposit.

Adsorption still requires gas reactants and vacant solid nodes, producing
surface adsorbates. Desorption produces vacant nodes and gas products; it now
accepts solid reactants in `gas_exposed` mode. Reactions can express other local
transformations. Diffusion conserves species **and phase** inventories, allowing
both adsorbate hops and solid-atom moves, including connected concerted moves.
Molecular formulas are not inferred from labels; chemical atom conservation is
still the input author's responsibility. One occupied node represents one
specified effective solid unit, not a resolved multi-atom molecule.

## Counting and rates

- `embeddings`: rate per ordered assignment of slot roles to distinct sites.
  Symmetric permutations count separately; choose this only when that rate
  convention is intended.
- `unique_site_sets`: rate per unordered set of matched sites. All matching
  slot assignments for that set are grouped; select an assignment uniformly
  when the event fires. This prevents a symmetric three-site pattern from
  automatically multiplying its rate by its permutations.
- `sites`: a single-site rule, counted once per eligible node.
- `undirected_pairs`: two-slot alias for unordered site sets.
- `reactive_sites`: two-slot historical convention; both endpoints can be
  centers. Count eligible centers, then choose a matching assignment at that
  center. This preserves the ZGB rate convention and is distinct from per-edge
  and per-embedding rates.

Rates remain constant or Arrhenius, with the declared counting unit. The direct
KMC clock uses the sum of **eligible group counts × process rate** before the
event. The rate convention is an explicit modeling choice, not inferred from
species or geometry.

## Thickness, outputs and limits

In `gas_exposed` mode, mean thickness is the volume-equivalent quantity
`solid_count × site_volume / lateral_area`, in angstrom. For generated crystals,
site volume is the crystal volume per candidate node, excluding added vacuum.
Custom inputs specify it explicitly. This includes the initial substrate and
is not the maximum z or a height-map average for a porous film. Growth rate is
the change in that mean thickness divided by elapsed time, in angstrom/s.
The residence-time height integral uses this same quantity.

`column_tops` continues reporting height in layers and growth rate in layers/s.
Summaries include `height_unit`, `growth_rate_unit`, `surface_sites`, geometry
metadata and the process counters. XYZ preserves actual cell periodicity,
chemical labels and site IDs. Output filenames remain `surface_*.xyz`,
`surface_summary.json` and `surface.log`.

Reaching the top candidate plane by solid incorporation reports `capacity_limit`,
not poisoning. Add candidate headroom before interpreting a longer growth run.
The general matcher rebuilds access and eligible events after each update to
correctly handle global changes in connectivity or sealed cavities. It is
intended first for modest systems and does not yet have the large-system
throughput of the specialized ZGB catalogs. `max_embeddings` (default 100000)
and `max_search_steps` (default 1000000), in model `parameters`, bound each full
catalog rebuild across all processes. Exceeding either limit raises an error;
the matcher never silently truncates events or alters their total rate.
