# General system configuration

`eda_analysis_compact_pairwise_linear_csv_only.py` retains its original defaults,
but `--system-config FILE.json` can describe another substrate/catalyst family
without editing Python source. The configuration uses JSON so NumPy remains the
only third-party dependency.

## Important scope

One invocation must contain structures with the same molecular connectivity and
atom partition. Run different substrate/catalyst families separately when their
connectivity differs, then compare their CSV outputs. Atom numbers below are
one-based XYZ indices and must be replaced with indices from the user's files.
Contiguous atom sets may be written compactly as a quoted inclusive range, for
example `"ligand": "2-62"`. Do not write bare `2-62`, which is invalid JSON.
Explicit arrays and mixed arrays such as `[1, "3-8", 12]` are also accepted.

The isolated catalyst reference must retain the catalyst atom ordering used at
the beginning of each TS. The isolated substrate may have another ordering; it
is graph-mapped to the configured substrate.

## Quinone + bis-quinoline example

```json
{
  "files": {
    "quinone_con_R": "quinone_con_R.xyz",
    "quinone_con_S": "quinone_con_S.xyz",
    "quinone_step_R": "quinone_step_R.xyz",
    "quinone_step_S": "quinone_step_S.xyz"
  },
  "reference": "quinone_step_R",
  "structure_metadata": {
    "quinone_con_R": {"mechanism": "concerted", "stereo": "R", "substrate_family": "quinone", "catalyst_family": "bis_quinoline"},
    "quinone_con_S": {"mechanism": "concerted", "stereo": "S", "substrate_family": "quinone", "catalyst_family": "bis_quinoline"},
    "quinone_step_R": {"mechanism": "stepwise", "stereo": "R", "substrate_family": "quinone", "catalyst_family": "bis_quinoline"},
    "quinone_step_S": {"mechanism": "stepwise", "stereo": "S", "substrate_family": "quinone", "catalyst_family": "bis_quinoline"}
  },
  "fragments": {
    "metal": [1],
    "ligand": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    "reactive_atoms": [13, 14],
    "substrate": [15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
  },
  "reaction": {
    "substrate_sites": [
      {"atom": 17, "label": "quinone_C1"},
      {"atom": 18, "label": "quinone_C2"}
    ]
  },
  "descriptors": {
    "substrate_plane": [15, 16, 17, 18, 19, 20],
    "ligand_planes": {
      "quinoline_left": [2, 3, 4, 5, 6],
      "quinoline_right": [7, 8, 9, 10, 11]
    },
    "local_fragment": {"graph_radius": 1},
    "catalyst_core": [1, 3, 7, 11]
  }
}
```

The short atom lists are illustrative only; use the full heavy-atom plane for
each quinoline when that is the intended chemical plane.

Run it with:

```bash
python eda_analysis_compact_pairwise_linear_csv_only.py \
  --system-config quinone_bis_quinoline.json \
  --cat-ref bis_quinoline_cat.xyz \
  --sub-ref quinone.xyz \
  --eda EDA.csv \
  --all
```

## Plot-ready strong relationships

The default strong-trend filter is structure-level `|r| >= 0.90`, `R2 >= 0.80`,
and `n >= 5`. Configure it with:

```bash
python eda_analysis_compact_pairwise_linear_csv_only.py \
  --system-config system.json --all \
  --linear-min-abs-r 0.85 --linear-min-r2 0.72 \
  --linear-min-n 5 --linear-max-results 30
```

Outputs include:

- `selected_linear_relationships.csv`: ranked fit summary and robustness flag;
- `selected_linear_plot_data.csv`: all selected points in long format;
- `linear_relationship_tables/*.csv`: one directly plottable table per fit.

`--linear-max-results 0` exports every relationship passing the thresholds.
Pairwise-difference correlations remain in the complete correlation output but
are deliberately excluded from the strong-trend plot tables because their rows
are not independent observations.
