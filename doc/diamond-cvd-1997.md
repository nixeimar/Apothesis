# Diamond (111) CVD: 1997 mechanism

The source is C. C. Battaile, D. J. Srolovitz and J. E. Butler, *A kinetic Monte
Carlo method for the atomic-scale simulation of chemical vapor deposition:
Application to diamond*, J. Appl. Phys. **82**, 6293-6300 (1997). The supplied
paper is archived at [the University of Michigan](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/70750/JAPIAU-82-12-6293-1.pdf).

The `diamond_cvd_1997` backend implements all twelve Table II reaction families,
including the CH3 and C2H2 pathways, hydrogen reactions, dimer reconstruction and
carbon incorporation. This requires a C/H covalent bond-state model, rather
than labeling an occupied node as generic solid. The generic surface matcher
and existing ZGB models remain separate and available.

**Status:** all Table II channels are present, and atom conservation, the
(111) three-carbon bridge and (100) dimer-insertion paths are tested. The paper's
Figures 3-5 have **not** been numerically reproduced or statistically validated.
Several atom-to-lattice embedding conventions are explicit implementation
choices; a complete kinetic table alone does not uniquely specify them.
Do not interpret a short local run as a reproduction of the published curves.

## Run the supplied case

Rebuild the executable, then run:

```sh
./apothesis --validate-input examples/diamond-cvd-1997.json
./apothesis --input examples/diamond-cvd-1997.json
```

The root `input.json` is a copy of this diamond case. In Qt Creator, select the
project directory as the working directory and use the desired input's absolute
path in `--input`. Output goes to the working directory.

The local case uses 8 surface carbon atoms, 40 starting substrate carbon atoms,
8 terminating H atoms, 288 candidate carbon positions, T_surface = 1200 K,
T_gas = 1200 K, and a 200000-event limit. It requests two seconds of physical
simulation time, but the event cap can stop it earlier. Always check
`stop_reason` and actual `time` in the summary.

Gas partial pressures are H = 0.04, H2 = 18, CH3 = 0.004 and C2H2 = 0.04 Torr,
as in the paper's Figure 3. Figure 4 instead uses C2H2 = 0.1 Torr.
The source does not specify a unique gas temperature for the molar-density
conversion; 1200 K is an explicit input assumption and can be changed.

## Explicit mechanism and units

There are 22 independent directional channels in `processes`:

| Table II row | Forward / reverse process |
| --- | --- |
| 1 | H abstraction from diamond / H2-mediated repassivation |
| 2 | H addition to diamond / H desorption |
| 3 | H-assisted CH2 removal / CH2 adsorption from CH3 with H release |
| 4 | CH3 chemisorption / CH3 desorption |
| 5 | C2H2 chemisorption / C2H2 desorption |
| 6 | H abstraction from C1 hydrocarbon / H2-mediated addition |
| 7 | H addition to C1 hydrocarbon / H desorption |
| 8 | H abstraction from C2 hydrocarbon / H2-mediated addition |
| 9 | H addition to C2 hydrocarbon / H desorption |
| 10 | Addition of CH3 to a C1 hydrocarbon; no reverse supplied |
| 11 | Conditional dimer opening / dimer closing |
| 12 | Ring completion and insertion; no reverse supplied |

The fictitious example rates in Table I are **not** used. Table II parameters
are transcribed in `scripts/create_diamond_1997_input.py` and written explicitly
to each JSON channel. Remove or zero a channel to remove its events. `event`
selects a covalent-state operation, `domain` selects diamond/C1/C2 chemistry,
and the rate data remain in the input. These operations are C/H-specific
mechanism kernels; they are not inferred from an arbitrary reaction string.

Forward and reverse coefficient laws follow the paper:

- k_forward = A T^n exp(-E/(R T)).
- k_reverse = C A T^n exp(-deltaS/R) exp((deltaH-E)/(R T)).

E and deltaH are in kcal/mol; deltaS is in cal/(mol K). The implementation uses
R = 1.98720425864083 cal/(mol K), converting to kcal where needed. Negative
enthalpies and entropies must retain their signs. Gas coefficients are converted
to event propensities using [gas] = (P_Torr/760)/(82.057366080960 T_gas), in
mol/cm^3. They are not treated as raw per-site Arrhenius rates.

The standard-state correction 1.016e-5 mol/cm^3 is applied to reverse gas-addition
rows 2, 4, 5, 7 and 9. Other rows use unity, including the purely surface dimer
conversion. This is the implementation's dimensional interpretation of the
paper's standard-state prescription, and is explicit in each process input.

H reactions are counted per eligible C-H bond or free carbon valence direction.
An adsorption coefficient is per available anchoring direction: feasible
molecular orientations share that rate uniformly instead of multiplying it.
Dimer transitions and ring closures are per eligible configuration; insertion
orientations share the configuration rate. These degeneracy conventions are
implementation assumptions and must be audited before claiming agreement with
published numerical rates.

The clock is the direct KMC pre-event clock, dt = -ln(u)/sum(rates).
All fast H and dimer channels remain explicit; no time filter is used.
A compensated long-double accumulator retains very small increments during
stiff dimer cycling instead of rejecting or dropping those events.

## Carbon geometry and hydrogen representation

The generator creates an actual diamond-cubic (111) lattice. Its orthorhombic
cell vectors, in conventional cubic coordinates, are a/2[1,-1,0],
a/2[1,1,-2], and a[1,1,1]. Twelve carbon sites form the cell. Each atomic plane
has 2 nx ny carbon sites; a full cell contains six planes. The substrate stops
on a plane with three carbon back-bonds and one outward dangling direction,
which is initially H-passivated. Two bottom atomic planes are fixed.

The default lattice constant is an explicit geometry parameter of 3.567 A;
thermal expansion is not modeled. Carbon coordinates remain fixed throughout
the simulation. Hydrogen atoms decorate carbon valence directions separately
from candidate carbon occupancy. Their displayed positions use a chosen C-H
length of 1.09 A along the associated lattice direction; this length is an
input assumption, not a bond length fitted from this paper.

Each hydrogen is a real exported H record with a carbon parent and a C-H bond;
it is not a floating Cu-like growth particle or an AB2 pseudo-unit. H records
use IDs after the candidate-carbon range, with four possible directional slots
per carbon. Carbon atoms are separately classified as diamond or chemisorbed
hydrocarbon. Only diamond carbon contributes to film growth rate.

The chemistry uses a sigma-bond valence representation. Internal hydrocarbon
pi bonds, nonbonded steric forces and atomistic relaxation are not evaluated.
A C2H2 embedding initially has one H on each of its two carbons. The same rigid
carbon lattice is used before and after incorporation; molecular orientations
are discrete. These choices permit the competing pathways but are not a
validated reconstruction of every atomistic intermediate in the publication.

## Bond formation, dimers and growth

Hydrocarbon groups are connected components of chemisorbed carbon sigma bonds.
CH3 and C2H2 can adsorb/desorb and exchange H. Row 10 adds a methyl to a C1 group.
No surface diffusion is introduced, consistent with the paper's approximation.

A radical carbon pair can form a new nearest-neighbor C-C bond under row 12.
The engine requires an existing five-bond path, so the added edge closes a
six-carbon ring. If the resulting connected hydrocarbon has at least two bonds
to the diamond film, its carbons are relabeled as diamond. This includes the
three-carbon bridge from one C2H2 and one CH3 on (111). Relabeling conserves all
C and H atoms and leaves their explicit bonds intact.

A dimer can close between two surface diamond atoms that share a vacant carbon
position, have suitable free valences and are each bonded to no more than two
carbon atoms before closure. Dimer opening requires a neighboring radical
hydrocarbon. Insertion into an opened dimer rotates/reassigns a C1 or C2 molecule
to candidate positions, conserving C and H. C1 insertion retains H on the
inserted carbon; C2 insertion retains the remaining CHy group, following the
row-12 product form. Every accepted configuration obeys the carbon valence cap.
Dimer bonds are stored explicitly and exported separately; their rigid-lattice
endpoint separation is not a relaxed physical dimer bond length.

Initial carbon is connected to the fixed substrate. Chemisorbed carbon is
attached by explicit C-C bonds, and every transition preserves that connection.
C atoms cannot exceed four sigma/H/dimer bond valences. The simulation does not
include twinning, sp2 formation or thermal atomic vibrations.

## Defects and reported observables

Carbon occupancy defines void connectivity from the top gas boundary. Hydrogen
decorations do not create an artificial impermeable layer of candidate carbon
sites. After incorporation, enclosed empty carbon positions bounded by stable
diamond are marked as trapped vacancies; cavities still bounded by reversible
hydrocarbon are not permanently frozen. H on a buried diamond parent can be
marked trapped when it lies in the grown region and no reversible hydrocarbon
boundary remains. Trapped states no longer react.

This is an explicit defect-classification convention, not yet calibrated
against Figure 4. The implementation's normalized defect fractions use
N_new_diamond + N_trapped_vacancies + N_trapped_H as the denominator. The raw
counts are also reported so other sampling-volume definitions can be applied.

Mean diamond thickness = N_diamond a^3/(8 lateral_area), in A. Average growth
rate is its increase from the starting substrate divided by elapsed time.
Multiplying A/s by 0.36 gives micrometers/hour. Hydrocarbon residence does not
count as incorporated diamond. The exact residence-time thickness average is
reported separately from the growth rate.

Files are `diamond_summary.json`, `diamond.log`, `diamond_initial.xyz`,
`diamond_final.xyz`, `diamond_trajectory.xyz` and `diamond_bonds.tsv`. XYZ contains
actual C/H element records, stable site IDs, fixed flags, chemical state, H
parent and trapped flag. The bond file distinguishes carbon sigma, C-H and
virtual dimer bonds. Geometry output can be disabled with `write_xyz: false`.

## Larger cases and validation

Generate 300 surface carbon atoms and enough candidate headroom for the paper's
roughly 60 grown planes:

```sh
python3 scripts/create_diamond_1997_input.py --output diamond-300.json \
  --nx 10 --ny 15 --cells 12 --temperature 1200 --gas-temperature 1200 \
  --acetylene-pressure-torr 0.04 --duration 10 --max-events 10000000
./apothesis --input diamond-300.json
```

This is a larger input, not a promise of reaching ten seconds or sixty planes
within that event cap. The dimer channels are extremely stiff and can dominate
step count; inspect achieved time and deposited atoms. Do not treat capped,
induction-dominated or capacity-limited runs as steady growth-rate measurements.
The implementation rebuilds its event list and is a correctness-first backend,
not yet the optimized original authors' code. An unlimited event count (zero)
can require substantial computation.

For Figure 4 use C2H2 = 0.1 Torr and four independent films at each temperature.
Figure 3 and Figure 4 do not use the same acetylene pressure. Meaningful
comparison also requires the paper's substrate area, adequate deposited
thickness, multiple seeds, a growth/defect sampling convention audit and tests
of orientation/multiplicity assumptions. Those ensemble validations remain
outstanding.

Tests cover the full kinetic-table transcription, initial (111) passivation,
forward/reverse rate units, adsorption degeneracy, explicit C/H inventories,
(111) six-ring three-carbon bridging, dimer closing/opening/insertion, stochastic
C/H conservation, substrate connectivity, valence, seed reproducibility, the
pre-event clock and compensated accumulation of stiff events.

## Completed local run

The supplied 1200 K case, seed 1234, reached its 200000-event cap at
0.2704825346 s. Two row-12 events incorporated six carbon atoms, increasing
diamond carbon from 40 to 46. Two carbon atoms remained chemisorbed and nine
hydrogen atoms were present. Dimer opening and closing each occurred 85636
times, illustrating the stiffness of the unfiltered mechanism. No trapped H or
vacancies were recorded in this small initial sample.

An independent audit of the final exported XYZ and bond file verified total C/H
material balances, C-C and C-H lengths under the periodic cell, valence limits,
and connection of every atom to the fixed substrate. The local rate estimate
is an initial/transient eight-site result, not a validated steady growth rate.

To render actual exported atoms and bonds (Matplotlib required):

```sh
python3 scripts/plot_diamond_cvd.py --results /path/to/diamond-run
```
