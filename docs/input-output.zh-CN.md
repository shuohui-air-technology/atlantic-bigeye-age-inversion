# 用自己的数据运行

[返回首页](../README.zh-CN.md) · [English](input-output.md)

按下表准备月度重量 CPUE 和逐次作业渔获与努力量这两份 CSV。每个有观测的月份填一行，按 `YYYY-MM` 从早到晚排列。例如可以输入 `2020-11`、`2020-12`、`2021-02`、`2021-04`。模型会推进中间没有观测的月份，不把这些月份当成零渔获。

```sh
.venv/bin/bet-age-inversion --monthly-csv /path/to/monthly.csv --operations-csv /path/to/operations.csv --output-dir /path/to/results
```

至少提供两个观测月，首末月之间至少跨两次月度转移，才能定义进入过程的两个端点。每个重量档的月度数值至少有一个正值。作业年月必须出现在提供的观测月中。如果某月有观测且确实是零渔获，保留该月并将四档明确填写为零。也可以加上 `--selectivity tropical` 或 `--selectivity southern`。

Windows 下将 `.venv/bin/bet-age-inversion` 换成 `.\.venv\Scripts\bet-age-inversion.exe`。

## 输入与输出

| 输入 | 字段与单位 |
| --- | --- |
| 各重量档的月度 CPUE | 每个观测月一行，四列按重量档排列，单位为**千克整鱼重（RW）／千钩**。仓库示例有 18 个观测月，其中没有 2017 年 9 月的输入行。 |
| 月度 CSV | `month`、`CPUE_RW_lt15GG`、`CPUE_RW_15to25GG`、`CPUE_RW_25to40GG`、`CPUE_RW_ge40GG`。后缀 `GG` 表示用去鳃去内脏后重量划档，CPUE 数值已经按整鱼重（RW）计算。 |
| 作业记录 | `Year`、`Month`、`Hooks`、`wBET15D`、`wBET15U`、`wBET25`、`wBET40`。年月须与已拟合月份对应，钩数须大于零，各列须等长。 |
| 商业重量档顺序 | 单尾 GG 重量低于 15、15 至不足 25、25 至不足 40、40 千克及以上。四个作业重量字段分别填写该次作业各档渔获的总 GG 重量，单位为千克。 |

作业记录按 ICCAT 的换算系数 1.09 将 GG 重量转为 RW。按作业记录计算月度输入时，先按月、按档汇总渔获整鱼重，再除以当月所有作业的总钩数并乘以 1,000。零渔获作业的钩数也计入分母。

输出按 0 至不足 3 岁、3 至不足 5 岁、5 岁及以上的顺序排列。 `group_counts` 是尾数当量，`group_cpue` 是每千钩的尾数当量，`group_proportions` 是 0 到 1 的比例。零渔获作业各组的尾数当量和 CPUE 都为零，比例未定义，在 Python 中为 `NaN`，在示例 JSON 中为 `null`。跨作业汇总时，先加总尾数当量，再计算组比例。汇总 CPUE 使用总钩数。计算观测核和 `expected_age`（平均年龄）时，17 岁及以上的终末组统一以 17 岁表示。

### CSV 格式示例

以下各展示一行随附的合成输入。要直接运行，请下载完整文件。

**月度重量 CPUE。** 月度文件必须使用以下五列和顺序。

```csv
month,CPUE_RW_lt15GG,CPUE_RW_15to25GG,CPUE_RW_25to40GG,CPUE_RW_ge40GG
2016-09,2.1098379145468718,1.6346179598108228,2.463877008379354,12.168338797683813
```

**作业渔获和钩数。** 作业文件须包含以下七列，也可包含其他列。

```csv
Year,Month,Hooks,wBET15D,wBET15U,wBET25,wBET40
2016,9,1000,0,0,10,80
```

完整文件：[合成月度 CSV](../examples/synthetic/synthetic_monthly_weights.csv)和[合成作业 CSV](../examples/synthetic/synthetic_operations.csv)。这两份文件可以直接传给上面的命令。

### 输出文件

| 文件 | 内容 |
| --- | --- |
| `fit_result.json` | 拟合系数、收敛信息，以及有渔获月份的年龄组比例均值 |
| `monthly_age_indicators.csv` | 各观测月的年龄组 CPUE 和比例。观测到零渔获的月份，比例留空。 |
| `operation_results.json` | 包含 `group_counts_fish_equivalents`、`group_cpue_fish_equivalents_per_1000_hooks` 和 `group_proportions`。数组按作业 CSV 的行序排列。 |

## Python 接口

三种选择性设置都用同一套 Python 接口拟合月度状态和分配作业渔获：

```python
from bet_age_inversion import allocate_operations, fit_monthly_states

# monthly_cpue 的单位是千克 RW／千钩，operations 中的重量单位是千克 GG。
fit = fit_monthly_states(
    monthly_cpue,
    months=monthly_month_labels,
    bin_columns=["wBET15D", "wBET15U", "wBET25", "wBET40"],
    selectivity="tropical",
)
age_groups = allocate_operations(fit["bundle"], fit["states"], operations)
```

数据读取和完整调用见[命令行实现](../src/bet_age_inversion/cli.py)。固定的大眼金枪鱼（BET）配置支持平坦选择性 `selectivity="flat"`，可选正态、对数正态或 Gamma 观测核。平坦选择性下的正态核还支持相邻重量档转移，`epsilon` 可取 `0.05` 或 `0.10`。热带和南部体长选择性情景使用截断到正体长的正态核，不启用档间转移。`sample="full"` 拟合全部输入月份。原来的 `F1`、`F2` 和 `F3` 划分只用于论文中的月份序列。

## 模型设置

月度时间网格、各折的月份划分、进入过程的端点和求解设置见 [model_configuration.json](../src/bet_age_inversion/model_configuration.json)。生长关系、体长离散程度、自然死亡率、重量档边界、GG/RW 换算和月龄网格见 [parameters.py](../src/bet_age_inversion/parameters.py)。选择性曲线参数见 [selectivity.py](../src/bet_age_inversion/selectivity.py)。不同观测时段沿用同一套大西洋 BET 生物学设定和重量档边界。进入过程在用户所给时段的首次和末次月度转移之间插值。求解器从零开始估计三个非负系数，依次为初始年龄状态的尺度、首个进入端点和末个进入端点。分配作业渔获时，使用所选设置下经体长选择性加权的入档概率和重量一阶矩（每尾鱼对该档重量的期望贡献）。
