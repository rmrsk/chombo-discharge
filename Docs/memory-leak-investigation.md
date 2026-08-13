<!--
SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research

SPDX-License-Identifier: GPL-3.0-or-later
-->

# Memory growth investigation — status and handoff

Working notes for the unbounded memory growth that killed a 512-core ItoKMC run with
SIGKILL (exit 9). **Not a finished analysis** — the bug is not yet fixed. This records what is
established, what was ruled out, and the reproducer, so the work can be picked up without
repeating it.

## Summary

An ItoKMC simulation grows its memory footprint without bound as it runs, at a **fixed grid and a
fixed particle count**. It is not a leak in the conventional sense: nothing is orphaned, so valgrind
and ASan/LSan report the process as clean. The memory is reachable, owned, and never released.

A restart from a checkpoint reproduces the identical simulation state using **10.4x less memory per
cell**, which is the sharpest single demonstration that the footprint is a function of run history
rather than of the state being represented.

## Minimal reproducer

`Exec/Tests/BrownianWalker/DriftDiffusion` — Ito particles drifting and diffusing. No chemistry, no
photons, no field solve, no KMC.

```bash
cd Exec/Tests/BrownianWalker/DriftDiffusion
make -j DIM=2 MPI=TRUE DEBUG=TRUE OPT=HIGH
mpirun -np 4 ./main2d.*.ex regression2d.inputs \
  Driver.max_steps=400 Driver.plot_interval=999999 Driver.checkpoint_interval=-1 \
  Driver.write_memory=true Driver.regrid_interval=999999
```

400 steps, roughly 90 seconds, grid frozen:

| step | rank cells | tracked (bytes) | RSS |
|------|-----------|-----------------|--------|
| 1    | 2,420     | 2,266,372       | --     |
| 121  | 2,420     | 2,266,372       | 45 MB  |
| 241  | 2,420     | 2,266,372       | 49 MB  |
| 361  | 2,420     | 2,266,372       | 51 MB  |
| 396  | 2,420     | 2,266,372       | 51 MB  |

**Chombo's tracked memory is byte-identical for all 400 steps. RSS grows 13%.** The growth is
therefore entirely in memory Chombo's tracker cannot see.

Caveat: the step report prints no particle count, so the population has not been verified constant.
`DriftDiffusion` should neither create nor destroy particles, but confirm this before treating the
+6 MB as pure leak.

## What is established

- **Real memory.** RSS grows with all instrumentation disabled, in windows where the particle count
  moves by less than 0.5%. Not an artifact of the memory tracker.
- **Particle-driven.** With the grid frozen, over ~350 steps:
  - `CdrPlasma/DeterministicAir` (no particles): +1.3%
  - `CdrPlasma/StochasticAir` (photons only, 363 -> 11,254 verified): **0.0%, byte-identical**
  - `ItoKMC/AirBasic` (Ito + photons): **+11%**
  - `BrownianWalker/DriftDiffusion` (Ito only): tracked flat, **RSS +13%**
- **Ito-specific.** Photon particles run the same `ParticleContainer` / `EBAMRParticleMesh` /
  `Copier` machinery and do not leak.
- **Untracked.** Chombo's counter covers `Arena` (BaseFab), Chombo `Vector`, `Pool` and BitSet. It
  does **not** cover `ParticleSoA` arenas (`std::aligned_alloc`) or `std::vector` buffers. In the
  minimal reproducer the tracked total does not move at all.
- **Enters through `advance()`**, not `regrid()`. `Driver::regrid` frees ~78 KB net per call.

## Ruled out (by measurement, not argument)

`regrid()` as the entry point; `CoarseInterpQuadCF`; `Vector<VolIndex>`; Chombo `Pool`; the Krylov
work vectors; `snapshotGuess` (matched by `freeWork`, -14,171,076 vs +14,363,172); a `numSends`
ratchet; `Copier` accumulation (live count flat at 651); MPI request-vector growth (bounded at
numProc-1); `TraceTimer` (inert — `CH_TIMER` unset prunes the root); the memory tracker itself;
ParmParse rereads; photon particles; `ParticleSoA` slack (~3 MB and stable).

## Remaining surface

Untracked allocations on the Ito path:

- `ParticleContainer::remap()` -> `gatherMoversToPool` / `distributeFromPool`
- `MoverPool` — `std::vector<std::map<PoolKey, ParticleSoA>>`, one map **per rank** (1024 of them at
  production scale, built and destroyed every remap)
- the remap exchange buffers `sendBuf` / `sflat` / `rflat` in `CD_ParticleContainerImplem.H`
- `ItoSolver` per-step deposition and interpolation

**Next step:** bracket **RSS** (read `VmRSS` from `/proc/self/status`) around `ParticleContainer::remap()`
in the minimal reproducer. Do **not** use `ReportUnfreedMemory` — it is provably blind here.

## Secondary findings (real, but not the cause)

1. **Copier exchange buffers grow without bound.** `BoxLayoutDataI.H:593/605` allocate send/receive
   buffers under `if (size > capacity)` with no shrink path — a per-`Copier` high-water mark held for
   the Copier's lifetime. Confirmed and quantified in AirBasic (~17% of RSS growth there), but absent
   from the minimal reproducer, so not the driver.
2. **`AMRMultiGrid::init()` runs on every solve** from `CD_EllipticSolverChainImplem.H:142` with its
   re-entry guard commented out (`AMRMultiGrid.H:1212`). `revert()` frees nothing despite the comment
   at `:139` claiming it does.
3. **`Vector<T>::operator=(Vector<T>&&)`** (`Vector.H:344`) never decrements the destination's count.
4. **`RefCountedPtr`** charges construction to `Derived`'s arena and destruction to `Base`'s. This is
   the source of the negative `MFHelmholtzEBBCFactory` byte counts seen in production reports.

(3) and (4) make the per-type arena numbers unreliable; treat them as indicative only.

## Harness

- **Command-line overrides.** `ChomboDischarge::initialize` passes `argv+2` to `ParmParse`, so stock
  inputs files can be used unmodified: `./main.ex positive2d.inputs Driver.max_steps=400 ...`
- **Freeze the grid** with `Driver.regrid_interval=999999`. Preferred over coarsening the mesh, which
  suppresses refinement entirely and produces a misleadingly flat result.
- **A tiny `dt` freezes the grid but also stops the physics** — with `max_dt=1E-16` no photons are
  produced, so a photon-solver control run that way tests nothing. Verify the particles exist.
- `Driver.write_memory=true` plus `Driver.write_unfreed_memory=true` gives a per-allocation-site
  breakdown in `pout.<rank>`; RSS prints in every step report.

## Methodological warning

Two failure modes cost most of this investigation:

1. **One-sided brackets.** A bracket that catches an allocation whose matching release lives in
   another function always reads as a leak. This produced three separate wrong conclusions
   (`regrid`, `intersectParticles`, `snapshotGuess`). Always locate the paired release.
2. **Unrepresentative samples.** Growth is bursty — typically 80% of calls allocate nothing. Dumping
   "the first call with a non-zero delta" repeatedly landed on trivial calls. Select on a *large*
   delta.

Differential controls (compare two configurations differing in one thing) were reliable throughout.
Mechanisms inferred from reading code were not.

## Instrumentation

The phase-bracket instrumentation is committed alongside this document, in `CD_Driver.cpp`,
`CD_ItoKMCGodunovStepperImplem.H`, `CD_EllipticSolverChainImplem.H`, `CD_FieldSolverGMG.cpp` and
`CD_McPhoto.cpp`. It is diagnostic scaffolding, gated on existing flags, and **should be reverted
before this branch is merged anywhere.**

The Chombo submodule was also instrumented (`AMRMultiGrid.H`, `Copier.H`, `Copier.cpp`,
`BoxLayoutDataI.H`). Those changes cannot be committed here — the submodule is a separate
repository — so they are saved as `Docs/memory-leak-chombo-instrumentation.patch`. **A checkout of
this branch will not have them**, and the submodule in the working tree where this was done is dirty.
