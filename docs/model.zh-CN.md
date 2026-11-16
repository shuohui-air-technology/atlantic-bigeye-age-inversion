# 模型公式

[返回首页](../README.zh-CN.md) · [English](model.md) · [参数和数据来源](../SOURCES.md)

**月度输入（论文式 1–2）。** 设 $Y^{\mathrm{GG}}_{i,b}$ 为第 $i$ 次作业在重量档 $b$ 的总 GG 重量，$H_i$ 为该次作业的钩数，$E_m$ 为月份 $m$ 的全部作业。先转为整鱼重，再求每千钩的重量 CPUE。

```math
Y^{\mathrm{RW}}_{i,b}=1.09\,Y^{\mathrm{GG}}_{i,b},
\qquad
I_{b,m}=1000\,
\frac{\sum_{i\in E_m}Y^{\mathrm{RW}}_{i,b}}{\sum_{i\in E_m}H_i}.
```

**观测核（论文式 5–6，补充方法 S1）。** 设 $f(\ell\mid a)$ 为年龄 $a$ 的体长密度，$w(\ell)$ 为体长对应的单尾整鱼重，$B_b$ 为档 $b$ 的整鱼重区间。体长选择性 $s(\ell)$ 同时进入入档概率和重量一阶矩的积分。平坦参照下 $s(\ell)=1$。

```math
p_b^s(a)=\int_{\ell:\,w(\ell)\in B_b}s(\ell)f(\ell\mid a)\,\mathrm{d}\ell,
\qquad
K_b^s(a)=\int_{\ell:\,w(\ell)\in B_b}w(\ell)s(\ell)f(\ell\mid a)\,\mathrm{d}\ell.
```

**月度拟合（论文式 7–8）。** 月度年龄状态按月推进并计入自然死亡和新进入个体。$n_j(t_m)$ 是月份 $m$ 第 $j$ 个年龄格的相对状态，最后一格为 17 岁及以上。预测值由各年龄格的状态与该档重量一阶矩相乘后加总。论文中的比例系数 $q$ 吸收到拟合系数中。

```math
\widehat I_{b,m}=\Delta a\sum_{j=0}^{204}n_j(t_m)\,K_b^s(a_j),
\qquad
\Delta a=\frac{1}{12}\,\mathrm{yr}.
```

**作业分配（论文式 9–10）。** 先将该次作业对应月份的年龄状态归一化为 $\pi_m$，再用各档的整鱼重量计算逐龄尾数当量 $c_i(a_j)$。

```math
\pi_m(a_j)=\frac{n_j(t_m)}{\sum_k n_k(t_m)},
\qquad
c_i(a_j)=\sum_b Y^{\mathrm{RW}}_{i,b}
\frac{\pi_m(a_j)\,p_b^s(a_j)}{\sum_k\pi_m(a_k)\,K_b^s(a_k)}.
```

将组 $g$ 内的尾数当量相加，得到组尾数当量 $N_{i,g}$，再按该次作业的钩数换算为每千钩的尾数当量。

```math
N_{i,g}=\sum_{a_j\in g}c_i(a_j),
\qquad
\mathrm{CPUE}_{i,g}=1000\,\frac{N_{i,g}}{H_i}.
```
