# ZGB 2023 validation run

Executed on 2026-11-04 with five independent trials at each y_AB, for both growth and catalytic reference: 90 runs total. Settings: 100 × 100 periodic lattice, initial height 10, reaction rate 100000 s^-1, standard ZGB attempt normalization, site-based reaction selection. Base seed 1234; the reproduction runner records every individual seed.

Growth cases were run to absorption (2000 s censoring limit); all 45 absorbed. Catalytic cases used a 100 s horizon. Figure 2 height is integrated exactly with KMC residence times until poisoning, following the requested definition.

## Table 2 comparison

| y_AB | New poisoning time ± SEM (s) | Paper (s) | Relative difference |
| --- | --- | --- | --- |
| 0.1 | 557.36 ± 30.41 | 529.9 | +5.2% |
| 0.2 | 168.67 ± 12.96 | 165.2 | +2.1% |
| 0.3 | 86.21 ± 3.73 | 104.2 | -17.3% |
| 0.4 | 46.50 ± 2.05 | 50.8 | -8.5% |
| 0.5 | 33.91 ± 1.39 | 33.9 | +0.0% |
| 0.6 | 24.97 ± 0.87 | 21.4 | +16.7% |
| 0.7 | 17.50 ± 1.27 | 17.5 | +0.0% |
| 0.8 | 12.83 ± 0.54 | 13.0 | -1.3% |
| 0.9 | 12.23 ± 0.37 | 11.7 | +4.5% |

## Interpretation

The growth model reaches mixed, height-blocked absorbing states across all sampled compositions, while the reference retains a reactive window. The new mean growth-rate curve has a broad maximum around y_AB = 0.5–0.6; for this seed set the largest sample mean is at 0.5. This is qualitative agreement, not exact agreement with the paper's maximum at 0.6.

Some poisoning times differ materially (notably y_AB = 0.3 and 0.6). The paper provides no per-trial datasets or error bars, so no invented statistical acceptance criterion is applied. Exact numerical reproduction of every figure remains unverified. The implementation records its pair counting, reaction-center selection, pre-event clock, seeds, censoring and averaging to make remaining differences inspectable.

Unit tests independently recompute event catalogs after every transition and verify reactant balance, one height increment per growth reaction, four-neighbor equal-height constraints, periodic wraparound, exact time integration, reaction selection, seed reproducibility and horizon handling. Tests passed with AddressSanitizer and UndefinedBehaviorSanitizer; executable integration tests passed.

## Regenerate

```sh
python3 scripts/reproduce_zgb_2023.py --binary ./apothesis --output /path/to/new-experiment
python3 scripts/plot_zgb_2023.py --results /path/to/new-experiment
```

See [model and reproduction guide](zgb-growth-2023.md). The full experiment retains the per-trial inputs, logs, final matrices, summary JSON files and a comparison CSV.
