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

An ItoKMC simulation grows its memory footprint without bound as it runs. It is not a leak in the
conventional sense: nothing is orphaned, so valgrind and ASan/LSan report the process as clean. The
memory is reachable, owned, and never released.

A restart from a checkpoint reproduces the identical simulation state using **10.4x less memory per
cell**, which is the sharpest single demonstration that the footprint is a function of run history
rather than of the state being represented.

**Resolved: the tracked growth was a phantom.** `overallMemoryUsage` ratcheted on every regrid
because Chombo's `Vector<T>::operator=(Vector<T>&& )` performs **no memory accounting at all** — it
never decrements the destination's element count. Every move-assign onto a non-empty `Vector`
therefore leaves the destination's old size on the books permanently, with no allocation behind it.
Fixing that one omission removes the entire per-regrid ratchet: growth per regrid goes from a steady
+47,664 bytes to **exactly zero**. See "The cause" below.

This does **not** explain the 512-core SIGKILL, which is an RSS event and is not driven by the
tracked counter. The real resident-memory ratchets are documented as secondary findings 1 and 5.

## Reproducer

`Exec/Examples/ItoKMC/AirBasic` with `example.inputs` and a developing discharge. The grid must be
allowed to regrid — that is the trigger (see below), and freezing it crashes a developing discharge.

```bash
cd Exec/Examples/ItoKMC/AirBasic
make -j DIM=2 MPI=TRUE DEBUG=TRUE OPT=HIGH
mpirun -np 4 ./main2d.*.ex example.inputs \
  Driver.max_steps=150 Driver.plot_interval=-1 Driver.checkpoint_interval=-1 \
  Driver.write_memory=true Driver.write_unfreed_memory=true AmrMesh.max_amr_depth=2
```

150 steps, roughly 5 minutes. Tracked memory grows **+655,620 bytes**, and it grows on **19 of 149
steps — every one of them a regrid step or its immediate neighbour. All other steps are
byte-identical.**

### `BrownianWalker/DriftDiffusion` is NOT a reproducer — retracted

The previous handoff proposed `DriftDiffusion` as the minimal reproducer on the strength of "tracked
flat, RSS +13%", flagging the unverified particle count as a caveat. **The caveat was the whole
story.** The step report does print the population (`BrownianWalkerStepper::printStepReport`), and it
is not constant:

| config | global particles, step 1 -> 400 | rank-0 local, step 1 -> 400 |
|--------|--------------------------------|------------------------------|
| `ppc=32` (default, superparticle merge/split on) | 18,691 -> 35,833 (**+92%**) | 330 -> 20,485 (**62x**) |
| `ppc=0` (merging disabled) | 91,222 -> 63,099 (**-31%**) | 324 -> 48,090 (**148x**) |

The blob drifts across the domain, so rank 0 starts nearly empty and fills up. The reported "13% RSS
growth" is rank 0 acquiring 20,000 particles, at a steady ~313 bytes per particle. Two controls
settle it:

- **Serial.** `-np 1`, `ppc=32`: the particle count nearly doubles (18,575 -> 35,572) and RSS moves
  **+192 kB in 400 steps** (0.2%). Flat.
- **Per rank.** `-np 4`: ranks 1, 2 and 3 are flat to within 0.4 MB across the run. Every kilobyte of
  the growth is on rank 0 — the one rank whose load grows.

`DriftDiffusion` also has **no field solve at all**, so it could never have exercised the mechanism
below. Its tracked memory is byte-identical for the same reason it is byte-identical in any run that
never regrids.

## The cause

`Chombo-3.3/lib/src/BaseTools/Vector.H:127`:

```cpp
inline
Vector<T>& operator=(Vector<T>&& a_invec)
{
  v=std::move(a_invec.v);          // no increment(), no decrement()
  return *this;
}
```

Compare the copy-assign directly beneath it, which brackets the assignment with
`decrement(size())` / `increment(size())`. The move overload has neither.

Walk the counter through `dst = std::move(src)`, with `dst` holding N elements and `src` holding M:

| | real memory | tracked counter |
|---|---|---|
| before | dst: N, src: M | N + M |
| after the move | dst: M, src: 0 | N + M (unchanged) |
| `src` destroyed | — | −0 (it is empty; `~Vector` decrements `size()`) |
| `dst` destroyed | M freed | −M |
| **net left on the books** | 0 | **+N, forever** |

So every move-assign onto a non-empty `Vector` inflates the reported unfreed memory by the
destination's prior size. Nothing leaks; the counter is simply wrong, and wrong in one direction
only, so it accumulates monotonically with the number of move-assignments performed — that is, with
run history.

`Vector<VolIndex>` dominates because the operator/stencil rebuild path move-assigns those vectors
heavily, once per regrid, at a size proportional to the irregular-cell count.

### Confirmation (differential control)

Add the missing `decrement(size())` before the move — saved as
`Docs/chombo-vector-move-tracking-fix.patch` — and rerun the same cases unchanged:

| case | tracked growth before | after |
|------|----------------------|-------|
| `CdrPlasma/DeterministicAir`, 50 steps | +518,942 B | +279,022 B |
| `ItoKMC/AirBasic`, 150 steps | +655,620 B | +209,636 B |
| per-regrid `Vector<VolIndex>` (steady state) | **+47,664 B, every regrid** | **0** |

The per-step picture in AirBasic is the clearest statement of the result. Before, 19 steps grew, all
of them regrid steps. After, there is a single one-off event at step 15 (+209,636 B — genuine, the
grid is still growing as the discharge develops) followed by a **perfectly balanced** alternation of
−21,024 / +21,024 that nets exactly zero. No monotone component remains.

The residual after the fix, measured on `DeterministicAir` over a steady-state window (steps 34–50,
8 regrids), is +15,684 B total — about 2 kB per regrid, concentrated in `Arena 12EBDataImplem`,
`Arena 13EBGraphImplem` and `Arena 16EBISLayoutImplem`, i.e. EBIS layout caching. AirBasic showed no
residual at all over the same kind of window, so this may just be the tail of grid development
rather than a ratchet. Unconfirmed either way, and three orders of magnitude below the original
symptom.

### What this means for the production failure

- The **10.4x restart discrepancy in unfreed memory is largely this artifact.** The phantom grows
  with the number of move-assignments the process has performed; a restart from a checkpoint has
  performed almost none, so it reports a far smaller figure for an identical physical state. That
  observation was the anchor of this whole investigation and it was substantially measuring a
  counter bug.
- The **SIGKILL is not explained by this**, and never could have been: the scheduler kills on RSS,
  and the tracked counter has no bearing on RSS. For real resident growth, see secondary findings 1
  (Copier exchange buffers) and 5 (`ParticleSoA` capacity), both of which are genuine high-water-mark
  ratchets in memory the tracker cannot see.

Note also that finding 3 below — the same `Vector` move-assign omission — was already recorded in
this document as a caveat that "makes the per-type arena numbers unreliable". It was not a caveat.
It was the entire tracked-memory symptom.

## What was established along the way

**The tracked memory ratchets once per regrid, and is otherwise perfectly flat.**

In the AirBasic reproducer above, the per-step net change in `overallMemoryUsage` is **exactly zero
on every step that does not regrid**, and positive on every step that does. Growth is not a slow
per-step drip; it is a discrete event tied to grid changes.

Differential control — hold everything fixed and vary only the regrid frequency:

| `Driver.regrid_interval` | regrids in 150 steps | tracked growth | per regrid |
|--------------------------|----------------------|----------------|------------|
| 2  | 12 | 655,620 B | ~54.6 kB |
| 4  | 9  | 512,964 B | ~57.0 kB |
| 8  | 7  | 452,900 B | ~64.7 kB |

Fewer regrids, less growth, at a roughly constant cost per regrid. (The scaling is not exactly
linear because a different regrid cadence also produces different grids.)

This accounts for the production symptoms:

- **The 10.4x restart discrepancy.** A restart from a checkpoint reproduces the state with no regrid
  history behind it, so none of the ratchet has accumulated. This is the single sharpest piece of
  evidence in the investigation and the regrid mechanism predicts it directly.
- **Fixed grid, fixed particle count, unbounded growth.** The *instantaneous* grid is fixed; the
  number of regrids performed is not. A 512-core run regridding every few steps for 10^5 steps
  performs ~10^4 regrids, at a per-regrid cost that scales with grid size (~55 kB here on a tiny 2-D
  depth-2 grid).
- **Why valgrind/ASan/LSan are clean.** Nothing is orphaned; the memory is reachable and owned.

It also resolves the contradiction in the previous handoff. "Enters through `advance()`, not
`regrid()`" is correct as far as it goes — the *allocation* happens during the field solve inside
`advance()` — but it is *triggered* by the preceding regrid. Bracketing `regrid()` alone therefore
showed nothing, which is why the trail went cold.

Note that `CdrPlasma/DeterministicAir` (+1.3%) and `StochasticAir` (byte-identical) were both run
with the grid frozen, so they were controls against the wrong variable.

## Ruled out (by measurement, not argument)

`CoarseInterpQuadCF`; `Vector<VolIndex>`; Chombo `Pool`; a `numSends` ratchet; `Copier` accumulation
(live count flat at 651); MPI request-vector growth (bounded at numProc-1); `TraceTimer` (inert —
`CH_TIMER` unset prunes the root); the memory tracker itself; ParmParse rereads; photon particles.

Also ruled out, and re-confirmed this session with a matched pair of brackets:

- **The Krylov work vectors / `snapshotGuess`.** On every regrid step `FREEWORK delta = -153,828`
  and the reallocation in `snapshotGuess` is `+153,828`. Exactly balanced, every time.
- **`ParticleSoA` arenas as the cause of the *tracked* growth.** They are allocated with
  `std::aligned_alloc` and are invisible to `overallMemoryUsage`. They do ratchet (see secondary
  finding 5) and they do inflate RSS, but they cannot move the tracked counter.

## Remaining surface

The tracked-memory question is closed. What is left is **real resident memory** — and the useful
sorting principle is *what a regrid does and does not reset*, because production regrids every few
steps.

**Reset by a regrid, therefore bounded by the regrid interval:**

- `ParticleSoA` capacity (secondary finding 5). `ParticleContainer::regrid` news fresh holders and
  destroys the old arenas — measured, the 408x over-allocation collapses to 1.0x.
- Copier exchange buffers (secondary finding 1), `BoxLayoutDataI.H:593/605`. Visible in the post-fix
  AirBasic trace as a ±21,024 alternation — down at the regrid, back up after — netting exactly zero.

Neither is a plausible cause of an out-of-memory kill on a run that regrids regularly. Both were
measured with the grid frozen, which overstates them badly.

**Not reset by a regrid — investigated and cleared:**

`EBISLevel::m_cache` (`Chombo-3.3/lib/src/EBTools/EBISLevel.{H,cpp}`) was the obvious remaining
suspect, because it is the one structure a regrid never touches. It is a
`std::map<DisjointBoxLayout, EBISLayout>`, one per EBIS level, owned by an `EBIndexSpace` that is
built once from the geometry and lives for the whole run. How it works:

- **Keys are address-ordered, not content-ordered.** `BoxLayout::operator<` ends in
  `return m_layout < rhs.m_layout`, and `m_layout` is a `RefCountedPtr<int>` with no `operator<`,
  only an implicit `operator T*()` — so the comparison is on raw pointers. Two structurally identical
  grids from different regrids are therefore different keys, and **a regrid can never hit this
  cache**, even when it reproduces the mesh exactly.
- **Entries are expensive.** `EBISLayoutImplem::define` builds a `LevelData<EBGraph>` with
  `nghost+1` ghosts and a `LevelData<EBData>` with `nghost` ghosts over the requested layout, both
  filled by `copyTo` from the level's master copies, plus a `LayoutData<EBISBox>` wrapping them. That
  is real geometry: connectivity, volume fractions, centroids, areas, normals.
- **Entries nest.** `EBISLayoutImplem` holds `Vector<EBISLayout> m_coarLevels`/`m_fineLevels`, and
  `setMaxCoarseningRatio` recursively calls `fillEBISLayout` on coarser levels, so one retained entry
  can pin entries in other levels' caches.
- **Eviction is refcount-driven.** `refreshCache()` erases entries at `refCount() == 1` (cache holds
  the only reference). It runs from `fillEBISLayout` whenever `m_cacheStale == 1`, which is set on a
  miss and cleared immediately — so the sweep runs on every miss and only on a miss. In
  chombo-discharge the holders are the `EBLevelGrid`s built in `PhaseRealm::defineEBLevelGrid`.

The failure mode to look for was a holder outliving its grids, which would keep an entry above
refcount one forever. **It does not happen.** Instrumented with per-level `entries` / `pinned` /
`hits` / `misses` (in the instrumentation patch, reported by `Driver::writeMemoryUsage`), over 250
steps of `CdrPlasma/DeterministicAir` with ~120 regrids:

| step | 20 | 40 | 80 | 120 | 160 | 200 | 240 |
|------|----|----|----|-----|-----|-----|-----|
| entries | 20 | 24 | 24 | 24 | 24 | 24 | 24 |
| pinned  | 16 | 18 | 18 | 18 | 18 | 18 | 18 |

Both saturate by step 40 and never move again. Eviction works, and the earlier +2 kB/regrid residual
was the tail of grid development, exactly as that caveat warned.

### The tracked-memory question is fully closed

With the accounting fix in place, over the same 250 steps only 41 of 249 step-deltas are non-zero,
and they decompose into:

- **Balanced oscillation** — ±12,624, ±11,296, ±11,584 pairs, the Copier exchange buffers falling at
  a regrid and rebuilding after. Net zero.
- **Three discrete jumps** — steps 31, 79 and 147, totalling ~3.9 MB of the 4.18 MB — each coinciding
  with the grid genuinely growing. The per-regrid `Regrid AmrMesh` allocation, a proxy for grid size,
  goes 255,628 -> 320,528 -> 490,424 bytes at exactly those regrids.
- **Nothing at all after step 169.**

What remains is a high-water mark in grid size — when the grid later shrinks (proxy back to 308,196)
the memory is not returned — not growth that is unbounded in time. There is no per-regrid component
left.

**Next step:** measure **RSS**, not `overallMemoryUsage`. The tracked counter undercounts by a wide
margin (9 MB tracked against 58 MB RSS on a 2-rank AirBasic run), and the SIGKILL is an RSS event.
Nothing in the tracked counter now points anywhere.

Worth deciding separately: whether to carry the `Vector` move-assign fix as a local Chombo patch, and
whether to offer it upstream. It is a two-line change and it makes the tracked counter usable again.

## Secondary findings (real, but not the cause)

1. **Copier exchange buffers grow without bound.** `BoxLayoutDataI.H:593/605` allocate send/receive
   buffers under `if (size > capacity)` with no shrink path — a per-`Copier` high-water mark held for
   the Copier's lifetime. Confirmed and quantified in AirBasic (~17% of RSS growth there), but absent
   from the minimal reproducer, so not the driver.
2. **`AMRMultiGrid::init()` runs on every solve** from `CD_EllipticSolverChainImplem.H:142` with its
   re-entry guard commented out (`AMRMultiGrid.H:1212`). `revert()` frees nothing despite the comment
   at `:139` claiming it does.
3. **`Vector<T>::operator=(Vector<T>&&)`** (`Vector.H:127`) never decrements the destination's count.
   **This turned out to be the cause of the entire tracked-memory symptom** — see "The cause" above.
   Left listed here because it was originally filed as a minor accounting caveat, which is precisely
   the misjudgement worth remembering.
4. **`RefCountedPtr`** charges construction to `Derived`'s arena and destruction to `Base`'s. This is
   the source of the negative `MFHelmholtzEBBCFactory` byte counts seen in production reports.

5. **`ParticleSoA` capacity is a per-patch high-water mark that is never released.**
   `ParticleSoA::reserve` returns early unless the request exceeds the current capacity, and
   `shrinkToFit()` — which exists and is correct — **has zero callers anywhere in the tree**. So a
   patch's arena is sized to the largest particle population it has ever held, for the lifetime of
   the run. Measured in `DriftDiffusion` (`-np 4`, `ppc=0`, 400 steps, grid frozen), summing the
   `Bulk` container over all patches:

   | | step 1 | step 400 |
   |---|--------|----------|
   | live particles | 91,222 | 63,099 (**-31%**) |
   | allocated capacity | 176,433 | 227,331 (**+29%**) |
   | capacity / live | 1.93 | **3.60**, still climbing |

   Per rank the ratchet is stark: rank 3 ends the run holding **123 live particles in capacity for
   50,133** — a 408x over-allocation, 5.7 MB for 14 kB of data — having peaked at 44,407 particles
   and given back nothing. Rank 3's RSS is identical at step 1 and step 400 despite losing 99.7% of
   its particles.

   **A regrid resets this completely**, which bounds it far more tightly than the table above
   suggests. `ParticleContainer::regrid` allocates fresh `LevelParticles` holders
   (`new LevelParticles(a_grids[lvl])`), destroying the old arenas outright. Measured on the same
   rank 3, same case, varying only the regrid interval:

   | | grid frozen | regrid every 10 steps |
   |---|---|---|
   | final capacity for ~120 live particles | 50,133 | **64** |
   | peak capacity/live ratio | 408x | 2.25x |

   So the ratchet is bounded by the *regrid interval*, not by the length of the run. The 408x figure
   comes from a configuration production never runs in — it was measured with the grid frozen, which
   the earlier methodology used as its standard control. At a production regrid cadence of a few
   steps this is close to a non-issue. Still worth a shrink policy on its own merits (with a
   hysteresis threshold such as `capacity > 2 * size`, rather than an exact `shrinkToFit` that would
   force a reallocation on the next append), but it is not a candidate for the SIGKILL.

   It is also invisible to `overallMemoryUsage` — these arenas use `std::aligned_alloc` — so it was
   never a candidate for the tracked growth either.

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

   There is a live instance of the inverse error still in the tree. In
   `CD_EllipticSolverChainImplem.H`, `snapshotGuess` opens its bracket at `snap0` *after* calling
   `allocateWork()`, so the bracket structurally cannot observe the allocation it was written to
   measure — it prints `0` on precisely the steps where the work vectors are reallocated. Check that
   a bracket encloses the call you think it encloses before believing a zero.
2. **Unrepresentative samples.** Growth is bursty — typically 80% of calls allocate nothing. Dumping
   "the first call with a non-zero delta" repeatedly landed on trivial calls. Select on a *large*
   delta. In the AirBasic reproducer the extreme case holds: 130 of 149 steps allocate *exactly zero*
   tracked bytes, and 100% of the growth is on the other 19.
3. **Unverified caveats are not minor.** The `DriftDiffusion` reproducer carried one explicit
   unverified assumption — that the particle count was constant — and that assumption was false by a
   factor of 62 on the measured rank. It invalidated the whole reproducer. Verify the caveat before
   building on the result, not after.
4. **Control against the right variable.** Freezing the grid was used throughout as the way to hold
   conditions fixed. Since regridding turns out to *be* the trigger, freezing it suppressed the very
   effect being hunted, and every "byte-identical" control run was measuring its absence.

Differential controls (compare two configurations differing in one thing) were reliable throughout.
Mechanisms inferred from reading code were not.

## Instrumentation

The phase-bracket instrumentation is committed alongside this document, in `CD_Driver.cpp`,
`CD_ItoKMCGodunovStepperImplem.H`, `CD_EllipticSolverChainImplem.H`, `CD_FieldSolverGMG.cpp` and
`CD_McPhoto.cpp`. It is diagnostic scaffolding, gated on existing flags, and **should be reverted
before this branch is merged anywhere.**

Added while establishing the above, and likewise to be reverted:

- `CD_Driver.cpp` — the step report printed RSS and unfreed memory as `std::ceil(bytes/MB)`. At 1 MB
  granularity the entire effect documented here is invisible: the AirBasic growth is 655 kB in 150
  steps, and it would have read as a flat `5(MB)` throughout. RSS now prints in kB and unfreed memory
  in bytes. **Any future measurement here needs byte or kB resolution** — the MB rounding is what
  made the tracked counter look "byte-identical" in cases where it was in fact moving.
- `CD_BrownianWalkerStepper.cpp` — `printStepReport` also sums `ParticleSoA` size and capacity across
  all patches (the `#soa` line), which is what secondary finding 5 is measured with.

The Chombo submodule was also instrumented (`AMRMultiGrid.H`, `Copier.H`, `Copier.cpp`,
`BoxLayoutDataI.H`). Those changes cannot be committed here — the submodule is a separate
repository — so they are saved as `Docs/memory-leak-chombo-instrumentation.patch`. **A checkout of
this branch will not have them**, and the submodule in the working tree where this was done is dirty.

The `Vector` move-assign fix is saved separately as `Docs/chombo-vector-move-tracking-fix.patch`,
because unlike the rest it is **not** scaffolding — it is a real bug fix and should not be reverted
along with the diagnostics. Apply both from `Submodules/Chombo-3.3`:

```bash
cd Submodules/Chombo-3.3
git apply ../../Docs/chombo-vector-move-tracking-fix.patch
git apply ../../Docs/memory-leak-chombo-instrumentation.patch   # diagnostics only
```
