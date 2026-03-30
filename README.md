# PDE MATLAB Experiments

这个仓库收集了一组用 MATLAB 实现的 Monge--Amp\`ere 方程数值实验原型，重点包括：

- 论文 Example 5.1 的 `vanishing moment + Legendre-Galerkin` 复现
- 多种从“较易辅助问题”到 Monge--Amp\`ere 的同伦原型
- 当前主线的 `u_0` 锚定同伦 Legendre-Galerkin 实验
- 算例汇总、批量对比、TeX/PDF 图集生成

大多数脚本都在物理区域 `(0,1)^2` 上工作，内部使用参考区间 `[-1,1]^2` 的 Legendre 基与 B-matrix 变换链。

## 目录概览

- `paper_example_5_1_legendre_demo.m`
  论文 Example 5.1 的基准复现脚本。
- `u0_anchor_monge_ampere_homotopy_demo.m`
  当前 `u_0` 锚定同伦主入口，也是这个仓库最近工作量最多的一条线。
- `u0_anchor_case_summary_demo.m`
  批量跑一组算例并输出总表。
- `u0_anchor_case_gallery_n_sweep_demo.m`
  为图集批量生成 paired `(N,p)` 的误差随 `N` 图表。
- `semilinear_to_monge_ampere_homotopy_demo.m`
  半线性到 Monge--Amp\`ere 的独立同伦原型。
- `poisson_to_monge_ampere_homotopy_bmatrix_demo.m`
  `Poisson -> Monge--Amp\`ere` 的实验脚本。
- `biharmonic_to_monge_ampere_homotopy_bmatrix_demo.m`
  `biharmonic -> Monge--Amp\`ere` 的实验脚本。
- `convex_t_biharmonic_homotopy_demo.m`
  与凸性控制相关的四阶同伦实验。
- `output/`
  所有数值结果、CSV、图片、TeX、PDF 默认都写在这里。

## 快速开始

### 1. 运行论文 Example 5.1 基准复现

```matlab
paper_example_5_1_legendre_demo
```

输出目录：

```text
output/paper_example_5_1/
```

### 2. 运行 `u_0` 锚定同伦单算例

```matlab
u0_anchor_monge_ampere_homotopy_demo('paper_example_5_1')
```

输出目录类似：

```text
output/u0_anchor_monge_ampere_homotopy_paper_example_5_1/
```

### 3. 运行 paired `(N,p)` sweep

当前仓库里常用的 paired 取法是：

- `N = 2,4,6,\ldots,32`
- `p = 3,5,7,\ldots,33`
- `t = 10^{-p}`

示例：

```matlab
opts = struct( ...
    'N_values', [2, 4:2:32], ...
    'p_values', [3, 5:2:33], ...
    'show_figures', false);
u0_anchor_monge_ampere_homotopy_demo('paper_example_5_1', opts)
```

这会生成：

- `u0_anchor_error_vs_N.png`
- `u0_anchor_N_convergence.csv`
- `u0_anchor_N_convergence.md`

## `u_0` 锚定同伦主线

主脚本：[`u0_anchor_monge_ampere_homotopy_demo.m`](u0_anchor_monge_ampere_homotopy_demo.m)

求解的是：

```math
t (u-u_0) + (1-t)\bigl(\det(D^2u)-f\bigr)=0,
\qquad t:1\to 0.
```

其中：

- `t=1` 对应锚定态 `u=u_0`
- `t=0` 回到目标 Monge--Amp\`ere 方程 `det(D^2u)=f`
- 默认边界条件始终与目标问题共享同一个 `g`

这条主线当前支持：

- 单算例 continuation
- paired `(N,p)` sweep
- 连续范数 `L^2 / H^1` 误差汇总
- 网格指标 `GridL2 / FDH1 / FDH2`
- Hessian 正定比例、Jacobian 条件数等诊断

## 这次更新了什么

最近一轮整理主要集中在 `u_0` 锚定同伦与图集工作流：

- 给 paired `(N,p)` sweep 加入了论文风格的连续范数列
  - `FinalContL2`
  - `FinalContH1`
- paired sweep 现在从 `(N,p)=(2,3)` 开始
- 批量图集 `u0_anchor_case_gallery` 已改成中文正文和中文图注
- 图集里每个算例都包含：
  - 方程 `u, f, g`
  - `u0_anchor_homotopy_error_curves.png`
  - `u0_anchor_error_vs_N.png`
  - 一张仿论文 Table 5.2 的结果表
- `quartic_convex` 目前使用了专门的分离型 anchor，以改善高 `N` 下的稳定性

## 批量脚本

### 算例总表

```matlab
u0_anchor_case_summary_demo
```

输出目录：

```text
output/u0_anchor_case_summary/
```

主要文件：

- `u0_anchor_case_summary.csv`
- `u0_anchor_case_summary.md`
- `u0_anchor_conditions_report.md`

### 图集与 PDF

```matlab
u0_anchor_case_gallery_n_sweep_demo
```

然后在 `output/u0_anchor_case_summary/` 下编译：

```text
u0_anchor_case_gallery.tex
u0_anchor_case_gallery.pdf
```

当前图集 PDF 路径：

```text
output/u0_anchor_case_summary/u0_anchor_case_gallery.pdf
```

## 输出目录约定

常见输出目录命名方式如下：

- 单算例：
  - `output/u0_anchor_monge_ampere_homotopy_<case>/`
- `N`-sweep：
  - `output/u0_anchor_monge_ampere_homotopy_<case>_N_sweep_<tag>/`
- 图集与总表：
  - `output/u0_anchor_case_summary/`

通常每次运行会产出这些文件中的一部分：

- `u0_anchor_history.csv`
- `u0_anchor_error_vs_t.png`
- `u0_anchor_homotopy_error_curves.png`
- `u0_anchor_error_vs_N.png`
- `u0_anchor_N_convergence.csv`
- `u0_anchor_N_convergence.md`

## 已知现象与说明

- `paper_example_5_1` 这类光滑、强凸、`f>0` 的算例通常非常稳。
- `quartic_convex` 虽然本身很光滑，但默认 Poisson anchor 对它不够贴合；当前高 `N` 成功主要依赖专门的分离型 anchor。
- `trig_convex` 属于更难的近退化/振荡型算例，当前主线还不能稳定说明“满足光滑条件就一定成功”。
- paired `(N,p)` 图不是纯固定 `t` 的空间收敛图；它把 `N` 增大和 `t=10^{-p}` 变小绑在一起看。

## 运行环境

建议环境：

- MATLAB
- Windows + PowerShell（当前仓库的批处理与调用记录主要在这个环境下完成）
- 若需要编译图集或文章 PDF，建议安装 `XeLaTeX`

## 推荐阅读顺序

如果第一次看这个仓库，建议按下面顺序：

1. `paper_example_5_1_legendre_demo.m`
2. `u0_anchor_monge_ampere_homotopy_demo.m`
3. `u0_anchor_case_summary_demo.m`
4. `u0_anchor_case_gallery_n_sweep_demo.m`
5. `semilinear_to_monge_ampere_homotopy_demo.m`

这样比较容易先建立“论文基准 -> 当前主线 -> 批量结果 -> 其他同伦原型”的整体脉络。
