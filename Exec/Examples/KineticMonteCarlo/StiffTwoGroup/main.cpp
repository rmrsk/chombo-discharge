// Minimal 0D driver: run many independent realizations of the chemistry in chemistry.json and
// print the mean total electron number at stop_time.
#include <limits>
#include "ParmParse.H"
#include <CD_ItoKMCJSON.H>

using namespace ChomboDischarge;
using namespace Physics::ItoKMC;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  ParmParse pp(argc - 2, argv + 2, NULL, argv[1]);

  int  numRuns;
  Real stopTime, maxDt;
  pp.get("num_runs", numRuns);
  pp.get("stop_time", stopTime);
  pp.get("max_dt", maxDt);

  Random::seed();

  auto      physics   = RefCountedPtr<ItoKMCJSON>(new ItoKMCJSON());
  const int numPlasma = physics->getNumPlasmaSpecies();
  const int numPhoton = physics->getNumPhotonSpecies();

  Vector<Real>     particles(numPlasma, 0.0), photons(numPhoton, 0.0), phi(numPlasma, 0.0);
  Vector<RealVect> gradPhi(numPlasma, RealVect::Zero);

  physics->defineKMC();

  Real sum = 0.0;
  for (int run = 0; run < numRuns; run++) {
    for (int i = 0; i < numPlasma; i++) {
      particles[i] = 0.0;
    }
    particles[0] = 1.0; // one electron in the reacting group ("e" is the first plasma species)

    for (Real t = 0.0; t < stopTime;) {
      const Real dt    = std::min(maxDt, stopTime - t);
      Real       newDt = std::numeric_limits<Real>::max();

      // Zero field: every rate in chemistry.json is a constant, so E is irrelevant.
      physics->advanceKMC(particles, photons, newDt, phi, gradPhi, dt, RealVect::Zero, RealVect::Zero, 1.0, 1.0);
      t += dt;
    }
    sum += particles[0] + particles[1]; // e + el
  }

  pout() << "mean total electron number = " << sum / numRuns << endl;

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
