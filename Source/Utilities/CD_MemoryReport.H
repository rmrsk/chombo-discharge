/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MemoryReport.H
 * @brief  Headers for utility functions that report memory usage
 * @author Robert Marskar
 */

#ifndef CD_MEMORYREPORT_H
#define CD_MEMORYREPORT_H

// Std includes
#include <iostream>

// Chombo includes
#include <REAL.H>
#include <Vector.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief Namespace containing utility functions for reporting memory usage in various formats
 */
namespace MemoryReport {

/**
 * @brief Print the global maximum and minimum memory usage. This includes peak memory and unfreed memory.
 */
void
getMaxMinMemoryUsage();

/**
 * @brief Print the global maximum and minimum memory usage. This includes peak memory and unfreed memory.
 * @param[out] a_maxPeak    Maximum peak memory (maximized over ranks)
 * @param[out] a_minPeak    Minimum peak memory (minimuzed over ranks)
 * @param[out] a_maxUnfreed Maximum unfreed memory (maximized over ranks)
 * @param[out] a_minUnfreed Minimum unfreed memory (minimized over ranks)
 */
void
getMaxMinMemoryUsage(Real& a_maxPeak, Real& a_minPeak, Real& a_maxUnfreed, Real& a_minUnfreed);

/**
 * @brief Get peak and unfreed memory usage per rank.
 * @details If using MPI, this will do an Allgather so that every rank can see usage.
 * @param[out] a_peak    Peak memory usage
 * @param[in]  a_unfreed Unfreed memory usage
 */
void
getMemoryUsage(Vector<Real>& a_peak, Vector<Real>& a_unfreed);

/**
 * @brief Print a per-allocation-site breakdown of the memory that is currently unfreed.
 * @details Each line of the report is a single allocation site, written as 'file:line: bytes', so
 * the output attributes the process' unfreed memory to the code that allocated it. This reports
 * only the calling rank's allocations, so it must be called on every rank whose usage is of
 * interest -- pout() is the natural sink because it is already per-rank.
 *
 * Only allocations that pass through Chombo's Arena are counted, i.e. BaseFab-backed mesh data such
 * as FArrayBox, EBCellFAB and LevelData. Particle arenas (ParticleSoA uses std::aligned_alloc) and
 * every std::vector-backed buffer are invisible to the tracker, so this is a breakdown of the mesh
 * memory only, not of the process' resident set.
 *
 * The report is a no-op unless the code was compiled with CH_USE_MEMORY_TRACKING.
 * @param[in,out] a_os Stream that the breakdown is written to.
 */
void
reportUnfreedMemory(std::ostream& a_os);

/**
 * @brief Get the calling process' current and peak resident set size.
 * @details This is the memory the operating system has actually given the process, read from
 * /proc/self/status (VmRSS and VmHWM). Unlike the rest of this namespace it counts every allocation
 * -- particle arenas, std::vector buffers, MPI internals and the binary itself -- not just the mesh
 * data that passes through Chombo's Arena. It is therefore the quantity that determines whether a
 * job is killed for exceeding its memory allocation, and the one to compare the tracked figures
 * against when deciding whether growth in those figures is real.
 *
 * Both outputs are set to zero on platforms without /proc, so a zero reading means "unavailable"
 * rather than "no memory in use".
 * @param[out] a_currentRSS Current resident set size (VmRSS), in bytes.
 * @param[out] a_peakRSS    Peak resident set size (VmHWM), in bytes.
 */
void
getResidentSetSize(Real& a_currentRSS, Real& a_peakRSS);

/**
 * @brief Get the global maximum and minimum resident set size over all ranks.
 * @details See getResidentSetSize() for what is being measured. The maximum is the number that
 * matters for out-of-memory failures, since the job dies when the worst rank exceeds its limit --
 * not when the average does.
 * @param[out] a_maxCurrent Maximum current RSS, in bytes (maximized over ranks).
 * @param[out] a_minCurrent Minimum current RSS, in bytes (minimized over ranks).
 * @param[out] a_maxPeak    Maximum peak RSS, in bytes (maximized over ranks).
 * @param[out] a_minPeak    Minimum peak RSS, in bytes (minimized over ranks).
 */
void
getMaxMinResidentSetSize(Real& a_maxCurrent, Real& a_minCurrent, Real& a_maxPeak, Real& a_minPeak);
} // namespace MemoryReport

#include <CD_NamespaceFooter.H>

#endif
