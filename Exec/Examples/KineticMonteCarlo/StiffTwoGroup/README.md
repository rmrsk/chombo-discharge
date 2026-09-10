## Exec/Examples/KineticMonteCarlo/StiffTwoGroup

Reproducer for [issue #728](https://github.com/chombo-discharge/chombo-discharge/issues/728): the
hybrid tau-leaping algorithms lose reactions when one reaction turns its reactant over faster than
the step.

The model is two electron groups, a slow branching reaction feeding a fast one that relaxes back:

```
e  -> el + el + M+     nu_i = 7.468e9  1/s
el -> e                nu_r = 1.741e11 1/s     (23x faster)
```

All rates are constants, so no transport data is needed and the field is irrelevant. Starting from
a single `e`, the total `e + el` grows at the dominant eigenvalue of `[-nu_i, nu_r; 2 nu_i, -nu_r]`
= `6.8987e9 1/s`. Every reaction is first order, so the moment hierarchy closes and the mean is
known exactly: **993.89** at `t = 1 ns`.

# Compilation

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 main```

# Running the example

```mpirun -np <num_proc> main2d.*ex example.inputs```

Each rank prints the mean total electron number over `num_runs` realizations to its `pout.*` file.
With 12 ranks and `num_runs = 20000` the sampling error on the mean is about 0.2%.

# What it shows

| configuration | mean at 1 ns | vs exact 993.89 |
|---|---|---|
| `ItoKMCJSON.algorithm=ssa max_dt=1e-9` | 994.93 | +0.10 % |
| as shipped (`hybrid_midpoint`, `max_dt=1e-11`) | 922.29 | **-7.20 %** |
| `max_dt=1e-12` | 995.08 | +0.12 % |
| `ItoKMCJSON.crit_num=500` | 336.36 | **-66.16 %** |
| `ItoKMCJSON.crit_num=500 ItoKMCJSON.SSA_lim=0.0` | 992.21 | -0.17 % |
| `ItoKMCJSON.crit_num=500 ItoKMCJSON.SSA_lim=1000.0` | 996.55 | +0.27 % |

SSA is correct, so the reference is not in doubt. Two independent problems are visible:

* **The step does not resolve the fast reaction.** At `max_dt = 1e-11` we have `nu_r * dt = 1.74`,
  so the relaxation is drawn `Poisson(133)` firings against a population of 76. The state is
  invalid, the step is rejected and halved, and the accepted step is conditioned on the random draw.
  Neither `KMCSolver::computeDt` nor `KMCSolver::getNonCriticalDt` catches it: both bound the *net*
  change per species, and for a fast intermediate in quasi-steady state production and loss nearly
  cancel.

* **The SSA/critical branch in `advanceHybrid` is chosen from a random draw.** `curDt` can be the
  exponential `dtCrit`, and `useSSA = (A * curDt < m_SSAlim)` then selects the branch from the
  realized value of that draw. The critical branch applies the pending critical firing; the SSA
  branch does not. The two pure regimes (`SSA_lim = 0` and `SSA_lim = 1000`) are both correct and
  the default mixture is not, which the last three rows show. `crit_num = 5` normally hides this
  because few reactions are ever critical.

See the issue for the full analysis.
