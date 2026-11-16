# Frozen synthetic example

[English README](../../README.md) · [中文版](../../README.zh-CN.md)

| File | Purpose |
| --- | --- |
| [synthetic_monthly_weights.csv](synthetic_monthly_weights.csv) | Input weight CPUE in kg round weight per 1,000 hooks for 18 observed months |
| [synthetic_operations.json](synthetic_operations.json) | Three invented operations with four bin weights in kg gilled-and-gutted weight and effort in hooks |
| [synthetic_operations.csv](synthetic_operations.csv) | The same three operations in CSV format for the installed command |
| [reference_result.json](reference_result.json) | Reference flat-selectivity fit summary |
| [reference_monthly_age_indicators.csv](reference_monthly_age_indicators.csv) | Reference monthly fish-equivalent CPUE and age-group proportions |

The monthly input comes from the fixed simulation with generating coefficients `[1.0, 0.8, 0.8]`, 5% relative noise and seed `20261017`. The example reads these values directly. It applies the fitted synthetic state to the three synthetic operations.

After installation, run the example from the repository root:

```sh
.venv/bin/python examples/run_example.py
.venv/bin/python examples/run_example.py --selectivity tropical --output-dir example_output/tropical
.venv/bin/python examples/run_example.py --selectivity southern --output-dir example_output/southern
```

The default run writes `fit_result.json`, `monthly_age_indicators.csv` and `operation_results.json` to `example_output/`. The first two files are checked against the bundled flat-selectivity references. The additional runs refit the same synthetic inputs under the tropical or southern selectivity curve. See [parameter sources](../../SOURCES.md) for their provenance.

这套示例的输入是 18 个月的合成月度数据和三条合成作业记录。默认运行会核对平坦选择性下的拟合结果和逐月指标，并用同一组合成年龄状态计算作业结果。热带和南部选项读取同一份合成输入重新拟合。
