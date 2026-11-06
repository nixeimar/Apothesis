# Explicit processes on a fixed 3D surface

Select `model.type: "surface_kmc_3d"`. This model supplies surface geometry and
an event engine; it does not supply chemical species or implicit adsorption or
reaction channels. Every event comes from `processes`. The
[explicit ZGB example](../examples/zgb-explicit-3d.json) uses
`surface_sites: "column_tops"` to preserve its original surface restrictions.
The [general geometry and matching guide](general-surface-geometry.md) describes
`gas_exposed` mode for other crystals, custom positions and multi-site patterns.

## Define adsorption

```json
{
  "id": "ab_adsorption",
  "type": "adsorption",
  "reactants": [
    {"species": "AB", "phase": "gas"},
    {"vacancy": true, "site": 0}
  ],
  "products": [{"species": "AB", "phase": "surface", "site": 0}],
  "rate": {"model": "constant", "parameters": {"value": 0.4}},
  "parameters": {"counting": "sites"}
}
```

`site` identifies a local event slot, not a global crystal-site ID. The engine
matches the pattern at eligible binding sites anywhere in the surface. Gas
terms specify a reservoir species and carry no site number. All species,
including reservoir species, must be declared in `species`. Use element `X`
for an effective chemical particle without a resolved atomic structure.

Dissociative adsorption uses two explicit vacant slots, producing an adsorbate
in each. In the example, B2 adsorption has constant rate 0.3 per eligible
**undirected** neighboring pair, equivalent to (1-y_AB)/2 with y_AB = 0.4.
There is no hidden `y_ab` parameter or neighbor normalization in this engine.
Changing composition requires changing the explicit AB and B2 rates together.

## Define growth

```json
{
  "id": "ab_b_growth",
  "type": "reaction",
  "reactants": [
    {"species": "AB", "phase": "surface", "site": 0},
    {"species": "B", "phase": "surface", "site": 1}
  ],
  "products": [
    {"species": "AB2", "phase": "solid", "site": 0},
    {"vacancy": true, "site": 1}
  ],
  "rate": {"model": "constant", "parameters": {"value": 100000}},
  "parameters": {
    "counting": "reactive_sites",
    "neighbors": "lateral_first",
    "same_height": true
  }
}
```

A `solid` product fills that actual 3D crystal site and exposes the empty node
above it for subsequent adsorption. It requires occupied first-neighbor support
beneath it. The neighboring slot becomes vacant at its original position.
To model a catalytic reaction instead, use a vacancy product at both slots,
and add any gas product as a separate gas term. To change the precursor or
product identity, change the species in the pattern; names are not hardcoded.
Removing a process removes its events completely. Independent process IDs keep
separate counters even when their patterns are identical.

In `column_tops` mode, every pattern must explicitly define the same slots on both sides: slot 0 for
single-site events, or slots 0 and 1 for pair events. Each slot appears once per
side. Write separate terms for pair slots rather than `coefficient: 2`.

## Rate counting

- `sites`: rate in s^-1 per matching single binding site.
- `undirected_pairs`: rate in s^-1 per eligible geometric edge, counted once.
  Either orientation may match. If both orientations match and lead to distinct
  outcomes, choose either with equal probability.
- `reactive_sites`: rate in s^-1 per center site that has at least one matching
  partner. Select a center uniformly, then a matching partner uniformly. Both
  ends can be centers, as in the historical ZGB implementation. This convention
  differs from a rate per pair and is an explicit part of the model definition.

The example reproduces the ZGB event eligibility and propensity conventions.
Catalog ordering differs between engines, so the same seed need not produce
the same trajectory. Scientific comparisons require repeated trials.

Constant rates use `rate.parameters.value`. Arrhenius rates use `prefactor`
(s^-1) and `activation_energy` (J/mol), with positive
`system.parameters.temperature` (K): k = prefactor exp(-E/(R T)). The resulting
k uses the process's declared counting convention. Environment-dependent
barriers are not implemented here. Model-specific rate expressions require a
new rate implementation; arbitrary parameters alone do not add physics.

## Desorption and diffusion

Desorption has surface reactants, vacant site products and at least one gas
product. It removes adsorbates; it does not remove the supporting solid.

Diffusion uses `type: "diffusion"`, reactants `surface species at slot 0` plus
`vacancy at slot 1`, and products `vacancy at slot 0` plus `same species at slot 1`.
Use `counting: "undirected_pairs"` and the same neighbor conditions. This moves
an adsorbate on the supported top surface. The separate atomic Cu backend still
supports solid-atom adsorption, diffusion and desorption.

## Scope and validation

The `column_tops` mode used in this example retains the paper's simple-cubic,
one-/two-site surface rules. `gas_exposed` mode supports generated SC/BCC/FCC/HCP/
diamond slabs and custom candidate coordinates, connected multi-site patterns,
configurable distance links, solid removal and simultaneous solid products.
See the [general geometry guide](general-surface-geometry.md) for counting,
connectivity, thickness units and performance limits. All modes use fixed
coordinates and effective chemical units. Molecular formulas and chemical atom
conservation are not inferred from species names.

Extra custom values in `parameters` are retained in the input document. They
affect dynamics only when the engine implements them. Structural misspellings,
undefined species, duplicate process IDs, duplicate slots and invalid rates are
rejected. Only constant and Arrhenius rate models are supported.

## Output

Files are written to the working directory, which is printed when a run finishes:

- `surface_initial.xyz`, `surface_final.xyz`, `surface_trajectory.xyz`: positions,
  stable site IDs, fixed flags, `chemical_species`, and role (`substrate`, `film`,
  `adsorbate`). These names distinguish generic runs from the old ZGB preset.
- `surface_summary.json`: per-process counts and evaluated rates, solid and
  adsorbate inventories, mean height, exact residence-time average height and
  average growth rate since the start.
- `surface.log`: time series with one counter column per process ID.

`output.write_xyz: false` disables XYZ files. `lattice_interval: 0` keeps initial
and final frames only. Filling the final available growth layer stops with
`capacity_limit`; it is not physical poisoning.

The old `cheimarios_zgb_2023_3d` and height-field presets remain available for
comparison. They retain their own input conventions and `zgb*` output filenames.
