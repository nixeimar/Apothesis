# Cheimarios 2023 ZGB growth reproduction

This predefined model implements Section 2 and Table 1 of:

N. Cheimarios, *Insights into the effect of growth on the Ziff–Gulari–Barshad model
and the film properties*, Modelling and Simulation in Materials Science and
Engineering **31**, 065007 (2023), DOI: [10.1088/1361-651X/ace276](https://doi.org/10.1088/1361-651X/ace276).

It is a separate height-field model. The fixed-site atomic Cu backend is not used
for this reproduction. No diffusion, standalone thermal desorption, arbitrary
chemical reactions or off-lattice relaxation are added to the 2023 mechanism.

## Run a case

After rebuilding:

```sh
./apothesis --validate-input examples/zgb-growth-2023.json
./apothesis --input examples/zgb-growth-2023.json
```

The supplied case uses a 100 × 100 square lattice, initial height 10, lateral
periodic boundaries, y_AB = 0.4 and reaction constant 100000 s^-1.
The default example stops at 100 s or earlier absorption. The root `input.json`
now defines this chemistry through [explicit 3D processes](explicit-surface-processes.md); this example selects
the original height-field reproduction.

For the nongrowing ZGB reference (including the Figure 3a composition):

```sh
./apothesis --input examples/zgb-catalytic-2023.json
```

## Surface state and rules

Each column stores an integer film height and exactly one adsorption state:
`vacant`, `AB`, or `B`. Deposited AB2 is represented by the column height,
independently of the adsorbate covering its top.

The four neighbors of a column are north, south, east and west, with periodic
wraparound. Diagonals are not neighbors.

1. **AB adsorption** occupies one vacant site; its height does not change.
2. **B2 dissociation** occupies two neighboring vacancies as B. In the growth
   variant, the two columns must have exactly the same height.
3. **Reaction** consumes neighboring AB and B, creates two vacant adsorption
   sites, and increases only the AB precursor column by one. In the growth variant,
   the reactants must have the same height before reacting. The new vacancies
   therefore occur at different heights.

`variant: "catalytic"` disables incorporation and the height restrictions, giving
the nongrowing reference. The reference does not acquire a film-height change.

## Rates and event counting

`model.parameters.y_ab` is in [0, 1]. `reaction_rate` defaults to 100000 s^-1.
The only three chemical events are predefined; user `processes` sections are not
accepted by this preset. Additional custom values can be stored in `parameters`.

A vacancy has AB adsorption propensity y_AB. Nearest-neighbor pairs are stored
**once per undirected edge**, not once per orientation. The implementation exposes
pair normalization explicitly because Table 1 does not describe that bookkeeping:

- `pair_convention: "zgb_attempts"` (default): an incident B2 molecule chooses a
  site and one of its four neighbors. Each eligible directed interaction has
  rate (1-y_AB)/4; combining both orientations gives (1-y_AB)/2 per undirected
  vacant pair. On an empty, flat lattice, total AB + B2 propensity is N.
- `pair_convention: "ordered_pairs"`: each eligible directed interaction has
  rate (1-y_AB), giving 2(1-y_AB) per undirected vacant pair. This intentionally
  changes the relative AB/B2 kinetics and is an audit option, not the default
  reproduction convention.

The standard normalization is also consistent with the
[July 2023 Apothesis example](https://github.com/nixeimar/Apothesis/blob/7d1d52b/src/input.kmc):
it uses AB = 0.4 and B2 = 0.15, together with a neighbor-count multiplier in the
adsorption implementation (0.15 = (1-0.4)/4).

`reaction_selection: "sites"` (default) follows the historical site-based
Apothesis selection: every occupied site having at least one eligible unlike
neighbor is active at `reaction_rate`; choose the center uniformly and then one
of its eligible neighbors uniformly. The same undirected reaction may be reached
from either endpoint. A center with several partners remains one active site.
This distinction matters when dissociation creates competing reactions.
`reaction_selection: "pairs"` is an explicit alternative: each undirected AB–B
pair is an event at `reaction_rate`. It is retained for auditing, not silently
substituted for site-based selection.

The KMC clock uses the **pre-event** total propensity. No event after the requested
end time is executed. The old post-event-rate clock in historical Apothesis is
not reproduced. Exact random-number streams and historical numerical bugs are
not needed or claimed for the new implementation.

Vacancies, equal-height vacant pairs, reaction bonds and reactive center sites
use indexed event catalogs. Only a changed column and its four neighbors need
refreshing. The engine does not scan the entire 100 × 100 lattice at every event.

## Observables and averaging

`zgb_summary.json` reports the seed, convention choices, stop reason and:

- `coverage_ab`, `coverage_b`, `coverage_vacant`: fractions of columns in each
  adsorption state at termination.
- `mean_height`: final spatial mean of all column heights, in ML.
- `time_average_height`: the residence-time-weighted spatial mean height over
  the elapsed run. It integrates the piecewise-constant KMC state exactly,
  rather than taking an unweighted average of output samples. This is the
  requested Figure 2 height observable: time average until poisoning.
- `mean_growth_rate`: net mean height increase divided by elapsed physical time,
  in ML/s; equivalently the time average of the deposition increments.
- `rms`: spatial standard deviation of final column heights.
- Event counters: AB adsorption, B2 adsorption, reaction and total events.

`zgb.log` contains sampled observables and counters. `zgb_final.csv` stores x, y,
height and top adsorbate, allowing Figure-3-style coverage and height maps.
`output.log_interval` defaults to 1 s; it affects saving, not event selection or
time integration. Samples are labeled with their actual event times.

An absorbing state has zero total propensity. For 0 < y_AB < 1 and positive
reaction rate, this means no vacancies and no allowed reactions. The recorded
time is the time of the final physical event, not an infinite wait. If the
requested time horizon or `max_events` is reached first, the run is censored and
must not be interpreted as a poisoning-time measurement. At y_AB = 0, dimer
adsorption can jam with isolated vacant sites; it need not yield exactly full
B coverage on a finite lattice.

## Repeat the published parameter sweep

```sh
python3 scripts/reproduce_zgb_2023.py \
  --binary ./apothesis \
  --output /path/to/new-zgb-experiment
```

Defaults:

- 100 × 100 columns, initial height 10, five independent trials per composition.
- y_AB = 0.1, 0.2, ..., 0.9 for both growth and catalytic reference.
- Growth cases continue until absorption, with a 2000 s censoring horizon.
  This extends the usual 100 s so slow cases such as y_AB = 0.1 can poison.
- Reference cases use a 100 s horizon, matching the paper's usual observation time.
- Deterministic, separately recorded trial seeds starting from 1234.

Choose additional points around 0.389 and 0.525 to resolve the reference transitions:

```sh
python3 scripts/reproduce_zgb_2023.py --binary ./apothesis \
  --output /path/to/finer-zgb-experiment \
  --variant catalytic --y 0.35,0.38,0.389,0.40,0.45,0.50,0.52,0.525,0.53,0.55
```

The runner saves every trial input and output, a manifest, and `comparison.csv`.
The CSV reports means, sample standard deviations and standard errors, censoring,
and comparisons against all nine poisoning times in the paper's Table 2.
Censored runs are not assigned a relative poisoning-time error.
A different empty experiment directory is required for each invocation.

## Plot the results

Plotting is optional and requires Matplotlib (`python3 -m pip install matplotlib`);
the solver and sweep runner have no Python plotting dependency.

```sh
python3 scripts/plot_zgb_2023.py --results /path/to/new-zgb-experiment
```

This generates `coverage.png`, `growth.png`, `poisoning-time.png` and (when y_AB = 0.2
was included) `surface-snapshots.png`. Error bars represent SEM across trials.
To add a catalytic snapshot at y_AB = 0.45 to the surface plot, run
`examples/zgb-catalytic-2023.json` in a separate directory and pass its CSV using
`--reference-surface /path/to/reference/zgb_final.csv`.

## Reproduction limits

The rules, input conventions and observables are explicit and testable. Agreement
with the published figures must be assessed from the saved data; implementing
the mechanism does not establish exact numerical reproduction. The paper does
not provide trial seeds, per-trial datasets or uncertainty estimates for Table 2.
There are also differences from the historical finite-rate reaction setup and
clock, documented above. The Table 2 comparison therefore reports deviations
without inventing a statistical pass/fail threshold or tuning rates to force a
match. Final coverage, final height and time-averaged height must not be conflated.

Tests check all four lateral interactions, exclusion of diagonals and unequal
heights, periodic wraparound, precursor-only height increments, adsorbate/material
balance, exact time integration, catalog consistency after every event,
reproducibility, pair normalization, reaction-center selection and time limits.
