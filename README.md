# Finite resources in black-hole growth

Core code accompanying **A Finite-Residence Bound on Obscured Black-Hole
Growth** (previous title: *Finite Resources Limit Obscuration and Capture in
Black-Hole Growth*). This repository contains the finite
residence/gas-budget calculation, a finite collisionless capture control, and
a relativistic, prescribed-input SIDM scattering experiment. It includes small
fixed inputs and scalar reference summaries; no manuscript sources, large data,
particle dumps, build products or caches are tracked.

## Quick start

Use Python 3.10–3.13 with NumPy, SciPy, Matplotlib and Numba. From a terminal:

```bash
git clone https://github.com/billkang-x/finite-resources-black-hole-growth.git
cd finite-resources-black-hole-growth
python -m pip install -r requirements.txt
python -m finite_resources.reproduce
python -m finite_resources.figures
python -m finite_resources.validation
```

Results go into `outputs/`. The first command verifies checksums, recomputes the
gas accounting and analytic incoming capture fraction, and checks the supplied
lightweight reference summaries. The second recreates the two Letter figures
and the larger inventory plot **from those saved summaries**. The third checks
elastic pair four-momentum conservation, frame transformations, analytic
Schwarzschild capture against adaptive geodesics, finite-layer cancellation,
and physical time units. It does not simulate a long-term halo.

Typical reference values are a fiducial gas cost of **2.565 × 10^9 solar masses**
at 50 Myr covering-weighted activity, a stationary retained fraction of
**8.996 × 10^-4** at dimensionless supply 2000 and half covering, compact
collisionless inventory **1.5185 × 10^-3 solar masses**, and a local capture
enhancement of **1.4305** from the two-seed summaries.

## Rerun the numerical experiments

A small local transport run (including first-run Numba compilation):

```bash
python -m finite_resources.inner_transport --case compact --cells 64 --dt 0.5 --duration 12000 --packets 6000 --sigma 1 --export-packets
```

For the Letter's baseline use the following command with each seed
`20261031` and `20261032`; set `--sigma 0` for its collisionless comparison:

```bash
python -m finite_resources.inner_transport --case compact --cells 128 --dt 0.25 --duration 96000 --packets 48000 --sigma 1 --seed 20261031 --export-packets
```

Space refinement: `--cells 256`; time refinement: `--dt 0.125`; particle
refinement: `--packets 96000`. The extended control uses `--case extended`.
The return-energy audit additionally uses 256 and 512 cells and the combinations
listed in `summary/inner_runs.json`. Stochastic streams can differ across Numba
or platform versions; assess conservation and ensemble differences rather than
requiring identical individual trajectories. The two-seed statistics are
diagnostics, not a confidence bound for arbitrary histories.

Reanalyze a newly written packet file by substituting its actual filename:

```bash
python -m finite_resources.packet_audit outputs/inner/compact_c128_dt0.25_T96000_N48000_x1_seed20261031_packets.npz
```

The geodesic-filtered finite orbit inventory is independently reproducible from
the small energy-distribution CSVs:

```bash
python -m finite_resources.orbit_control --case compact
python -m finite_resources.orbit_control --case extended
```

Defaults use the reported 12/16/384 energy/angular/radial quadrature orders.
Try `--energy-order 4 --angular-order 4 --radial-order 96` for a faster control.
The capturable mass on an orbit is initial inward flux times its radial period.
Uniform initial phase gives depletion `min(t/P, 1)`, without restoring particles
after absorption. The scalar JSON includes an independent spatial phase-volume
integral and mass balance. Binary quadrature outputs are generated locally.

To refit a nonnegative energy distribution rather than using the supplied
640-bin result:

```bash
python -m finite_resources.phase_space --case compact
python -m finite_resources.inner_transport --case compact --distribution outputs/preparation/compact_distribution.csv
```

This optional bounded least-squares fit can be much more expensive than the
quick checks. It fits density and radial pressure together, reports residuals,
and exits with failure if the declared moment tolerances are exceeded. Use
`--output-dir` on any CLI to select another output folder. All commands should
be run from this repository's root.

## What the calculations establish

- `bounds.py`: a mass element has finite column-time yield for outward,
  nonreturning trajectories with declared opacity, radius and radial-speed
  bounds. Time structure and angular structure are unrestricted within that
  class. The sampling bound also requires the stated sampling-measure cap.
- `gas_budget.py`: named supercritical exterior examples, a conservative
  threshold corner, continuity versus printed density normalization, matched
  finite-layer cancellation, and stationary fuel/wind retention accounting.
- `orbit_control.py`: finite collisionless capture inventory in a frozen
  potential, matched to the Schwarzschild barrier at 100 GM/c^2.
- `inner_transport.py`: local elastic scattering, geodesics and separate input,
  capture, return and resident ledgers. The duration 96000 GM/c^3 is approximately
  **13.1 hours**, for a **10^5 solar-mass** seed.

The 100 Myr incoming-energy comparison keeps the measured input fixed. It is a
conditional energy ceiling, not an evolved 100 Myr halo or a proof that all
SIDM growth models fail. The preparations use a frozen gas monopole, reflecting
outer wall and declared unresolved inner density continuation. Gas and dark
matter have separate budgets. See [PROVENANCE.md](PROVENANCE.md) and the paper
for the precise assumptions and scope.

## Inputs, outputs and citation

`inputs/` contains only two final density slices, their nonnegative energy
distributions and physical parameters. `summary/` contains sampled curves and
scalar/block summaries used for comparison. Their checksums are recorded in
`INPUT_SHA256.json`. Full formation histories, original particle-level outputs
and manuscript TeX/PDF files are intentionally absent. The core simulations
generate their own packet records if requested; all generated files remain
local and are ignored by Git.

## Revised Letter: CEERS proxy and residence-history analysis

The v1.1.0 additions reproduce the revised Letter's public-catalogue application,
conditional branch scan and residence-scaling figure. The original capture and
gas modules remain available. A quick run needs only the original dependencies
and the shipped **aggregate** count inputs:

```bash
python -m finite_resources.revision
python -m unittest discover -s tests -v
```

Outputs go to `outputs/revision/`: exact count intervals, a uniform-total-incidence
Beta posterior, conditional resource caps, the analytically profiled two-origin
mixture likelihood, labeled-history controls, and `residence_scaling.pdf/png`.
The figure uses full-precision calculated rates; annotations have opaque
backgrounds and leaders so plotted lines do not obscure their equations.
This quick command does **not** rerun selection on the external survey catalogue.

To reconstruct the selection from all 87,345 released CEERS DR1 rows, install
FITS support and download `ceers_cat_v1.0.fits.gz` from the
[CEERS DR1 catalogue directory](https://web.corral.tacc.utexas.edu/ceersdata/DR1/Catalog/).
Keep it outside the repository or inside the ignored `data/` directory:

```bash
python -m pip install -r requirements-observation.txt
python -m finite_resources.ceers_selection /path/to/ceers_cat_v1.0.fits.gz
python -m finite_resources.joint_constraints --selection-json outputs/ceers/public_selection_reconstruction.json
```

The default loader checks the paper-input SHA256. It reconstructs the four-band
detection parent separately from the colour-selected numerator and also reports
five-band, redshift-window and mass-window controls. Generated object records
remain local. The reference result is **N=1242, K=48, F=3.8647343%**, with an
exact 95% count interval of **[2.8630452%, 5.0916453%]**. The nominal Fiducial
disjoint-sector cap is **3.1146239%**; it lies inside that interval. Its central
value fails the declared necessary-condition comparison; this is **not a 95%
statistical rejection**, an SIDM-origin measurement or an official completeness
correction. Both the general maximum-cost bound and the more restrictive
disjoint fuel-plus-column application are output with their assumptions.

Additional interfaces:

```bash
python -m finite_resources.residence_history
python -m finite_resources.residence_history --input /path/to/labeled_history.json
python -m finite_resources.trajectory_qeff /path/to/wind_state.npz --threshold-km-s 1
python -m finite_resources.residence_scaling
```

`residence_history` integrates every visit of each unique parcel, counts external
material once, and compares short, outward, returning and fresh-supply histories.
`trajectory_qeff` extracts the column and radial streamline kernel from **the same
outward cells** of a user-supplied axisymmetric snapshot. A single snapshot is a
streamline diagnostic; it does not supply a time-resolved mass-weighted R_Q.
No snapshot or catalogue is bundled. Input schemas, provenance, statistical scope
and the distinction between the two budget formulas are documented in
[docs/REVISED_LETTER.md](docs/REVISED_LETTER.md).

Please cite the accompanying Letter when available. Background: Dong-Biao Kang,
*Seed-mass dependence of baryon-assisted black-hole growth in self-interacting
dark matter haloes*, Astrophysical Journal (accepted),
[DOI:10.3847/1538-4357/aeb5e1](https://doi.org/10.3847/1538-4357/aeb5e1).
That DOI belongs to the accepted article, not to this software or the Letter.
Software attribution is provided in `CITATION.cff`. Code is released under MIT.

## 中文简要说明

下载仓库并安装 `requirements.txt` 后，先运行
`python -m finite_resources.reproduce` 核查主要数值，再运行
`python -m finite_resources.figures` 生成图；结果保存在 `outputs/`。
运行 `python -m finite_resources.validation` 检验核心物理算法。
上述快速命令使用小型输入与汇总表，不会重新运行完整粒子实验；如需重新
计算，请使用上面的 `inner_transport` 与 `orbit_control` 命令。

最新修订新增内容可运行 `python -m finite_resources.revision`，包括 CEERS
代理计数、条件性分支扫描、历史核对照和驻留增强图。完整公开目录的选样
重建另用 `ceers_selection` 命令；需自行下载公开 FITS 文件并安装
`requirements-observation.txt`。这里的公开颜色代理分析未采用官方完备率，
也不从颜色计数推断 SIDM 起源比例。大型数据、轨迹快照、缓存和论文文件
均不上传。
