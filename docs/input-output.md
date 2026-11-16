# Use your own data

[Back to README](../README.md) · [简体中文](input-output.zh-CN.md)

Prepare a monthly weight CPUE CSV and an operation catch and effort CSV using the fields below. Give each observed month one row, sorted in increasing `YYYY-MM` order. For example, `2020-11`, `2020-12`, `2021-02` and `2021-04` are a valid calendar. The model advances through the missing months without treating them as zero-catch observations.

```sh
.venv/bin/bet-age-inversion --monthly-csv /path/to/monthly.csv --operations-csv /path/to/operations.csv --output-dir /path/to/results
```

Supply at least two observed months, with the first and last separated by at least two monthly transitions so that both entry endpoints can be defined. Each bin needs at least one positive monthly value. Operations must belong to the supplied observation months. If a month was observed with zero catch, keep that row and explicitly use four zeros. The alternative selectivity flags also work with the two CSV inputs.

On Windows, replace `.venv/bin/bet-age-inversion` with `.\.venv\Scripts\bet-age-inversion.exe`.

## Inputs and outputs

| Input | Fields and units |
| --- | --- |
| Monthly bin CPUE | One row per observed month and four columns in **kg round weight (RW) per 1,000 hooks**. The bundled example has 18 observations and no row for September 2017. |
| Monthly CSV | `month`, `CPUE_RW_lt15GG`, `CPUE_RW_15to25GG`, `CPUE_RW_25to40GG`, `CPUE_RW_ge40GG`. The suffix `GG` describes the bin boundaries. The CPUE values are already in RW. |
| Operation records | `Year`, `Month`, `Hooks`, `wBET15D`, `wBET15U`, `wBET25`, `wBET40`. Dates must match fitted months, the number of hooks must be positive, and columns must have equal length. |
| Commercial bin order | <15, 15–<25, 25–<40, ≥40 kg gilled-and-gutted (GG) per fish. The four **operation fields** are total **kg GG** in each bin. |

For operation records, the model converts GG to RW using the ICCAT factor 1.09. To calculate monthly input from operation records, sum each bin's catch in kg RW for the month, divide by the month's total hooks, including zero-catch operations, and multiply by 1,000.

The age groups are 0–<3, 3–<5 and ≥5 years, in that order. `group_counts` contains fish equivalents, `group_cpue` contains fish equivalents per 1,000 hooks, and `group_proportions` contains fractions from 0 to 1. Zero-catch operations have zero fish equivalents and CPUE and undefined proportions (`NaN` in Python, `null` in the example JSON). Across operations, sum fish equivalents before computing age proportions and use total hooks when computing pooled CPUE. The 17+ state is represented by age 17 in the kernel and in the `expected_age` summary.

### CSV examples

Each excerpt shows one row from the synthetic inputs. Use the linked full files to run the command.

**Monthly weight CPUE.** The monthly file requires these five columns in this order.

```csv
month,CPUE_RW_lt15GG,CPUE_RW_15to25GG,CPUE_RW_25to40GG,CPUE_RW_ge40GG
2016-09,2.1098379145468718,1.6346179598108228,2.463877008379354,12.168338797683813
```

**Operation catch and effort.** The operation file requires these seven columns and can contain additional columns.

```csv
Year,Month,Hooks,wBET15D,wBET15U,wBET25,wBET40
2016,9,1000,0,0,10,80
```

Full files: [synthetic monthly CSV](../examples/synthetic/synthetic_monthly_weights.csv) and [synthetic operation CSV](../examples/synthetic/synthetic_operations.csv). You can pass those files to the command above.

### Output files

| File | Contents |
| --- | --- |
| `fit_result.json` | Fitted coefficients, convergence information, and mean age-group shares over months with positive catch |
| `monthly_age_indicators.csv` | Age-group CPUE and fractions by observed month. Shares are blank for months with zero catch. |
| `operation_results.json` | Contains `group_counts_fish_equivalents`, `group_cpue_fish_equivalents_per_1000_hooks`, and `group_proportions`. Arrays follow the row order of the operation CSV. |

## Python API

The public Python entry points use the same monthly fit and operation allocation for all three selectivity settings:

```python
from bet_age_inversion import allocate_operations, fit_monthly_states

# monthly_cpue contains kg RW per 1,000 hooks; operations contains kg GG.
fit = fit_monthly_states(
    monthly_cpue,
    months=monthly_month_labels,
    bin_columns=["wBET15D", "wBET15U", "wBET25", "wBET40"],
    selectivity="tropical",
)
age_groups = allocate_operations(fit["bundle"], fit["states"], operations)
```

See the [command-line implementation](../src/bet_age_inversion/cli.py) for data loading and the complete call. The fixed BET configuration supports `selectivity="flat"`, with normal, lognormal or gamma kernels. The flat normal kernel also supports adjacent-bin transfer at `epsilon=0.05` or `0.10`. The tropical and southern length-selectivity scenarios use the positive-truncated normal kernel without bin transfer. `sample="full"` fits all supplied months. The original `F1`, `F2`, and `F3` folds apply only to the paper calendar.

## Model settings

The monthly grid, folds, entry knots and solver settings are recorded in [model_configuration.json](../src/bet_age_inversion/model_configuration.json). Growth, length dispersion, natural mortality, bin boundaries, GG/RW conversion and the monthly age grid are implemented in [parameters.py](../src/bet_age_inversion/parameters.py). The logistic parameters are in [selectivity.py](../src/bet_age_inversion/selectivity.py). The same Atlantic BET biological settings and bin definitions apply to all input calendars. Entry is interpolated between the first and last calendar transitions of the supplied period. The solver starts from zero and estimates three nonnegative coefficients in the order initial state scale, first entry endpoint, last entry endpoint. The operation calculation uses length-weighted bin-entry probabilities and first weight moments for the selected scenario.
