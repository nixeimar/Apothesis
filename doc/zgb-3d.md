# ZGB chemistry on explicit 3D lattice particles

This backend advances a real `Atomistic::Surface`: a volume of fixed crystal
coordinates, first-neighbor links and occupied/empty solid sites. Adsorbate states
are stored separately for each 3D site. It does not invoke the height-field
`ZGB::Simulation` to evolve the state or merely extrude a height-map trajectory.

The top-surface restrictions of the 2023 paper are deliberately retained. Thus it
is a volumetric fixed-site realization of that mechanism, not an off-lattice or
fully resolved molecular model.

## Run locally

The preset example selects this backend; the root `input.json` now uses the
[generic engine with explicit processes](explicit-surface-processes.md):

```sh
./apothesis --validate-input examples/zgb-growth-3d.json
./apothesis --input examples/zgb-growth-3d.json
```

The predefined preset is [`examples/zgb-growth-3d.json`](../examples/zgb-growth-3d.json).
It creates a 20 × 20 × 30 simple-cubic candidate lattice (12000 sites), with ten
initially occupied layers (4000 substrate particles), one fixed bottom layer,
and twenty empty growth layers. The reaction rate is 100000 s^-1 and y_AB = 0.4.
The lattice constant is a schematic coordinate scale of 1 angstrom; it is not a
calibrated copper crystal or chemical bond length.

A paper-sized run uses `repetitions: [100, 100, 30]`, retaining ten initial solid
layers and twenty growth layers. Set `output.write_xyz: false` for large ensembles
if full-volume snapshots are unnecessary.

## Actual occupancy changes

An exposed binding site is an empty 3D node directly above the highest solid
particle in its vertical first-neighbor chain. It has a real occupied neighbor
beneath it. Buried cavities and unsupported initial overhangs are rejected in
this paper-specific model; the separate general geometry builder still exists.

- **AB adsorption:** records `AB` at a vacant binding node. It does not create a
  solid atom and does not change the underlying solid occupancy.
- **B2 adsorption:** records `B` at two vacant, same-z lateral first-neighbor
  binding nodes. The four lateral interactions are derived from the 3D graph;
  diagonal and vertical pairs are excluded.
- **Reaction:** neighboring AB and B at the same z are consumed. In growth mode,
  the AB node becomes an occupied solid site labeled `AB2`. A new binding node
  becomes active immediately above this deposited particle. The B node stays
  empty and becomes available for adsorption at its original height.

`variant: "catalytic"` performs the adsorption/reaction chemistry on a flat 3D
substrate without adding solid particles. A rough catalytic reference input is
rejected, rather than silently changing the reference model's interactions.

Every deposited solid unit has an occupied first-neighbor support beneath it.
Each adsorbate also has a solid first neighbor. Initial contiguous substrate
chains and fixed bottom layers keep the film connected; no floating particles
or detached growth clusters are introduced by this mechanism.

## Chemical coarse graining

The paper defines AB, B and AB2 as chemical states but does not supply a molecular
basis, molecular orientations, bond lengths or coordinates for constituent atoms.
Therefore one point represents one **effective growth unit or adsorbate**, not
one constituent atom of a multi-atom molecule. No fictitious A/B coordinates are
invented. This distinction is important when counting exported particles.

The predefined model accepts the pseudo-element `X` for a generic substrate.
XYZ uses `X` for generic film/adsorbate particles, together with explicit
`chemical_state` and `role` properties. The original atomic Cu backend and
chemical-element validation remain available separately.

A reaction replaces two adsorbate particles with one solid AB2 unit. Particle
count decreases by one, but symbolic atom inventory is conserved:

- A inventory = number of AB adsorbates + deposited AB2 units.
- B inventory = AB adsorbates + B adsorbates + twice deposited AB2 units.

Resolving each molecule into constituent atoms would require an additional
molecular structure model and physical chemical identities; this backend does
not claim to provide that information.

## Rates and observables

Rate conventions, site-based reaction selection and the correct pre-event KMC
clock match the documented [2023 preset](zgb-growth-2023.md). Only the underlying
state representation differs. Optional `pair_convention` and `reaction_selection`
audit settings are retained and recorded in the summary.

Mean film height is derived from **occupied 3D solid sites / lateral site count**.
It is not an integer height array driving events. The time-averaged height is
integrated with exact residence times up to poisoning or the specified limit.
Coverages count adsorbates per lateral surface chain.

There is a finite number of empty crystal layers. If growth fills the highest
candidate site, the run stops with `capacity_limit`, `capacity_limited: true`
and `poisoned: false`. This is a geometry-capacity cutoff, not physical poisoning.
Increase the number of z repetitions and rerun. The default headroom is ample
for typical short 2023 examples but is not an unlimited-growth promise.

## Files

- `zgb3d_summary.json`: representation, geometry, rates, seed, stopping reason,
  solid inventory, coverages, mean height, time average and event counts.
- `zgb3d.log`: time series of the same principal observables.
- `zgb3d_initial.xyz`, `zgb3d_final.xyz`: explicit initial and final particles.
- `zgb3d_trajectory.xyz`: initial, intermediate and final extended XYZ frames.

XYZ contains Cartesian positions, the periodic cell, stable 3D site IDs, fixed
flags, `chemical_state` and `role`. States are `substrate`, `AB2_solid`,
`AB_adsorbate` and `B_adsorbate`. A vacancy is not emitted as a particle.
Color particles by `chemical_state` to distinguish deposited solid from surface
adsorbates. The trajectory is not a complete restart of the candidate graph/RNG.

`output.lattice_interval` controls intermediate frames in seconds. Zero saves
only initial/final trajectory frames. `output.write_xyz: false` disables geometry
files; `log_interval` still controls numerical output.

## Repeat trials using the 3D backend

```sh
python3 scripts/reproduce_zgb_2023.py \
  --binary ./apothesis --representation 3d --growth-layers 20 \
  --output /path/to/new-zgb-3d-experiment
```

This selects `cheimarios_zgb_2023_3d`, builds a 3D slab for each trial and reads
its 3D summaries. It does not run the 2D engine and render the results afterward.
The runner disables XYZ output for these potentially large ensembles. Capacity
cutoffs are counted separately in `capacity_limited_trials`; such trials must
not be interpreted as poisoning-time measurements. The existing plotting script
accepts the numeric comparison CSV from either representation.

## Validation

Tests independently project the volume **only as an oracle** to compare eligible
four-neighbor event sets against the height-field paper model after each event.
Production 3D dynamics never use that projection. Tests also check creation of
an actual solid site, movement of the binding site to a higher coordinate,
separate chemical state/solid occupancy, first-neighbor support, fixed substrate
protection, symbolic A/B inventory, time integration, seed reproducibility,
periodic pairs, diagonal exclusions and finite growth capacity.

Matching eligibility rules does not establish exact reproduction of every
published numerical figure. The shared model's documented statistical and
historical-clock limitations still apply. Raw 3D trial outputs should be compared
with the paper separately from the earlier height-field experiments.

### Completed 3D ensemble

The 100 × 100 × 30 explicit-site ensemble was run with five independent trials
at each y_AB = 0.1, 0.2, …, 0.9 for both growth and catalytic variants (90 runs).
All 45 growth trials reached an absorbing state, with no capacity cutoffs. The
catalytic trials used a 100 s observation horizon; 21 of 45 reached absorption
within that horizon. Remaining trials are censored, not poisoned.

For growth at y_AB = 0.4, the five-run mean poisoning time was 59.93 s
(sample standard deviation 11.57 s), versus 50.8 s in the paper's Table 2.
The mean residence-time-averaged height was 10.9685 layers and the final spatial
mean was 11.0600 layers. Across the nine compositions, mean poisoning-time
differences from Table 2 ranged from −12.2% to +18.0%; this is not a claim of
exact numerical reproduction.

The default 20 × 20 local example deposited 430 AB2 units in 1457 events,
ending at 28.7970 s with 4430 total solid units. An independent check of its
exported coordinates found solid first-neighbor support for every solid unit
and adsorbate. The final mean height was 11.0750 layers and its residence-time
average was 10.8884 layers.

To render actual exported coordinates (Matplotlib required):

```sh
python3 scripts/plot_zgb_3d.py --results /path/to/run-directory
```
