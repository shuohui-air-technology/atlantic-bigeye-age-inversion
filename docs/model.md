# Model equations

[Back to README](../README.md) · [简体中文](model.zh-CN.md) · [Parameter and data sources](../SOURCES.md)

**Monthly input (paper Eqs. 1–2).** Let $Y^{\mathrm{GG}}_{i,b}$ be the catch weight of operation $i$ in bin $b$, $H_i$ its hook count, and $E_m$ all operations in month $m$. Convert catch weight to round weight and calculate weight CPUE per 1,000 hooks.

```math
Y^{\mathrm{RW}}_{i,b}=1.09\,Y^{\mathrm{GG}}_{i,b},
\qquad
I_{b,m}=1000\,
\frac{\sum_{i\in E_m}Y^{\mathrm{RW}}_{i,b}}{\sum_{i\in E_m}H_i}.
```

**Observation kernel (paper Eqs. 5–6 and Supplementary Method S1).** Let $f(\ell\mid a)$ be the length density at age $a$, $w(\ell)$ the round weight of one fish at length $\ell$, and $B_b$ the bin's round-weight range. Length selectivity $s(\ell)$ weights both the bin-entry probability and the first weight moment. The flat reference uses $s(\ell)=1$.

```math
p_b^s(a)=\int_{\ell:\,w(\ell)\in B_b}s(\ell)f(\ell\mid a)\,\mathrm{d}\ell,
\qquad
K_b^s(a)=\int_{\ell:\,w(\ell)\in B_b}w(\ell)s(\ell)f(\ell\mid a)\,\mathrm{d}\ell.
```

**Monthly fit (paper Eqs. 7–8).** The age state advances monthly with survival and new entry. Here $n_j(t_m)$ is the relative state for age cell $j$ in month $m$, with ages 17 and above in the last cell. Predicted weight CPUE sums the state times the bin's first weight moment across ages. The fitted coefficients absorb the proportionality factor $q$ in the paper.

```math
\widehat I_{b,m}=\Delta a\sum_{j=0}^{204}n_j(t_m)\,K_b^s(a_j),
\qquad
\Delta a=\frac{1}{12}\,\mathrm{yr}.
```

**Operation allocation (paper Eqs. 9–10).** Normalize the month's age state to obtain $\pi_m$. Use each bin's round weight to compute fish equivalents at age, $c_i(a_j)$.

```math
\pi_m(a_j)=\frac{n_j(t_m)}{\sum_k n_k(t_m)},
\qquad
c_i(a_j)=\sum_b Y^{\mathrm{RW}}_{i,b}
\frac{\pi_m(a_j)\,p_b^s(a_j)}{\sum_k\pi_m(a_k)\,K_b^s(a_k)}.
```

Sum fish equivalents within group $g$ to obtain $N_{i,g}$, then use the operation's hook count to obtain fish equivalents per 1,000 hooks.

```math
N_{i,g}=\sum_{a_j\in g}c_i(a_j),
\qquad
\mathrm{CPUE}_{i,g}=1000\,\frac{N_{i,g}}{H_i}.
```
