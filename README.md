# Atlantic bigeye tuna age inversion

[简体中文](README.zh-CN.md) · [Input and output guide](docs/input-output.md) · [Model equations](docs/model.md) · [Sources](SOURCES.md)

This Python package estimates age-group catch indicators for Atlantic bigeye tuna (*Thunnus obesus*) from four commercial weight bins. It fits a monthly age state to weight-based catch per unit effort (CPUE) and uses that state to allocate each operation’s catch to ages 0–<3, 3–<5, and ≥5 years. Outputs include fish equivalents, age-group proportions, and fish equivalents per 1,000 hooks.

## Quick start

Run the bundled synthetic example from the repository root with Python 3.12.

```sh
python -m venv .venv
.venv/bin/python -m pip install -c constraints.txt .
.venv/bin/python examples/run_example.py
```

For Windows PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -c constraints.txt .
.\.venv\Scripts\python.exe examples/run_example.py
```

The example prints this line after it fits the flat-selectivity model, checks the reference values, and allocates catch by age.

```text
flat: fitted 18 synthetic months and allocated 3 synthetic operations.
```

The example saves its results in `example_output/`. The package requires Python 3.10 or later. We verified the synthetic example with Python 3.12 and the pinned versions of NumPy and SciPy. If those versions are unavailable for your Python version, run `python -m pip install .` in your virtual environment. See [the synthetic example](examples/synthetic/README.md) for the input sources and other selectivity options.

## Use your own data

Use two CSV files. The first gives monthly CPUE by weight bin (**kg round weight per 1,000 hooks**). The second gives each operation’s total catch in each bin (**kg gilled-and-gutted weight**) and its hook count. Supply the months you observed. The age state advances through months without observations.

```sh
.venv/bin/bet-age-inversion --monthly-csv /path/to/monthly.csv --operations-csv /path/to/operations.csv --output-dir /path/to/results
```

The [input and output guide](docs/input-output.md) lists the required columns, date formats, output fields, and Python API. The software reads and writes files locally.

## From bin weights to age-group CPUE

The monthly fit links the age state to catch weight through the weight-bin observation kernel. For operation $i$, let $c_i(a_j)$ be the fish equivalents at age $a_j$. Sum those equivalents within group $g$, then divide by the hook count $H_i$ and multiply by 1,000 to obtain fish equivalents per 1,000 hooks (paper Eq. 10).

```math
N_{i,g}=\sum_{a_j\in g}c_i(a_j),
\qquad
\mathrm{CPUE}_{i,g}=1000\,\frac{N_{i,g}}{H_i}.
```

See the [model equations](docs/model.md) for the weight conversion, selectivity weighting, monthly fit, and allocation by operation. The guide gives the equation numbers from the paper and links to the parameter sources.

## Repository layout

| Path | Contents |
| --- | --- |
| [src/bet_age_inversion/](src/bet_age_inversion/) | Installable Python package, model configuration, and command-line entry point |
| [examples/](examples/) | Synthetic runner, frozen inputs, and reference outputs |
| [docs/](docs/) | Input and output guide, Python API, and equations in English and Chinese |

## Data access and citation

The repository includes synthetic example data. For third-party public datasets used in the study, see the links to the original sources in [Sources](SOURCES.md).

The repository does not include the study’s commercial fishing records, vessel summaries, or coefficients fitted to the commercial data. For access to the study data, contact the corresponding author, [Peng Lian (lianpeng@cafs.ac.cn)](mailto:lianpeng@cafs.ac.cn).

Use [CITATION.cff](CITATION.cff) to cite the software. Record the version or commit you used and cite its Zenodo version DOI when available. Also cite the companion manuscript *Inferring Coarse Age Structure from Commercial Weight Bins: An Exploratory Approach for Atlantic Bigeye Tuna*, by Peng Lian, Haoran Yu, Jingcheng Li, Junzhou Shang, Lei Xu, Xiaomei Wang and Wenbo Yang.

## License

The software is available under the [MIT License](LICENSE).
