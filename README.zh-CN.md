# 大西洋大眼金枪鱼年龄反演

[English](README.md) · [输入与输出](docs/input-output.zh-CN.md) · [模型公式](docs/model.zh-CN.md) · [参数和数据来源](SOURCES.md)

这个 Python 包用四个商业重量档的数据，估算大西洋大眼金枪鱼（*Thunnus obesus*）各年龄组的渔获指标。模型利用月度重量 CPUE（单位捕捞努力量渔获量）拟合年龄状态，再将每次作业的渔获分配到 0 至不足 3 岁、3 至不足 5 岁、5 岁及以上三组。输出包括尾数当量、各组比例和每千钩尾数当量。

## 快速开始

在仓库根目录用 Python 3.12 运行合成示例。

```sh
python -m venv .venv
.venv/bin/python -m pip install -c constraints.txt .
.venv/bin/python examples/run_example.py
```

Windows PowerShell 使用以下命令：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -c constraints.txt .
.\.venv\Scripts\python.exe examples/run_example.py
```

完成平坦选择性拟合、参考值核对和作业分配后，程序会显示以下消息。

```text
flat: fitted 18 synthetic months and allocated 3 synthetic operations.
```

结果保存在 `example_output/`。包要求 Python 3.10 或更高版本。合成示例已在 Python 3.12 下用随附的 NumPy 和 SciPy 版本验证。如果所用 Python 版本无法安装这些依赖版本，可在虚拟环境中运行 `python -m pip install .`。输入来源和其他选择性选项见[合成示例说明](examples/synthetic/README.md)。

## 用自己的数据运行

运行时需要两份 CSV。第一份是各重量档的月度 CPUE（**千克整鱼重／千钩**）。第二份记录每次作业各档的总渔获重量（**千克去鳃去内脏重**）和钩数。按实际观测月份提供数据，年龄状态会继续推进到没有观测的月份。

```sh
.venv/bin/bet-age-inversion --monthly-csv /path/to/monthly.csv --operations-csv /path/to/operations.csv --output-dir /path/to/results
```

必需列名、日期格式、输出字段和 Python 接口见[输入与输出说明](docs/input-output.zh-CN.md)。文件均在本地读写。

## 从重量档到年龄组 CPUE

月度拟合通过重量档观测核，将年龄状态与各档渔获重量联系起来。设 $c_i(a_j)$ 为第 $i$ 次作业在年龄 $a_j$ 上的尾数当量。将年龄组 $g$ 内的尾数当量相加，除以该次作业的钩数 $H_i$，再换算成每千钩尾数当量（论文式 10）。

```math
N_{i,g}=\sum_{a_j\in g}c_i(a_j),
\qquad
\mathrm{CPUE}_{i,g}=1000\,\frac{N_{i,g}}{H_i}.
```

重量换算、选择性加权、月度拟合和作业分配的公式见[模型公式](docs/model.zh-CN.md)。该页标出了论文式号和参数来源。

## 仓库目录

| 目录 | 内容 |
| --- | --- |
| [src/bet_age_inversion/](src/bet_age_inversion/) | 可安装的 Python 源码、模型配置和命令行入口 |
| [examples/](examples/) | 合成示例脚本、固定输入和参考结果 |
| [docs/](docs/) | 中英文输入输出说明、Python 接口和公式 |

## 数据获取与引用

仓库随附合成示例数据。论文使用的第三方公开数据只提供原始来源链接，见[参数和数据来源](SOURCES.md)。

仓库不包含研究使用的商业作业记录、船舶汇总数据或基于商业数据拟合的系数。如需获取研究数据，请联系通讯作者 [Peng Lian（lianpeng@cafs.ac.cn）](mailto:lianpeng@cafs.ac.cn)。

软件引用信息见 [CITATION.cff](CITATION.cff)。请记录所用版本或提交编号。如果该版本已有 Zenodo DOI，请引用该版本 DOI。配套论文另行引用，题名为 *Inferring Coarse Age Structure from Commercial Weight Bins: An Exploratory Approach for Atlantic Bigeye Tuna*，作者为 Peng Lian、Haoran Yu、Jingcheng Li、Junzhou Shang、Lei Xu、Xiaomei Wang 和 Wenbo Yang。

## 许可证

本软件采用 [MIT 许可证](LICENSE)。
