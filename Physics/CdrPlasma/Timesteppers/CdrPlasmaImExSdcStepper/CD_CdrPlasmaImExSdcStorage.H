/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaImExSdcStorage.H
  @brief  Declaration of scratch storage for CdrPlasmaImExSdcStepper
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaImExSdcStorage_H
#define CD_CdrPlasmaImExSdcStorage_H

// Our includes
#include <CD_CdrPlasmaImExSdcStepper.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    class CdrPlasmaImExSdcStepper::CdrStorage
    {
    public:
      CdrStorage();
      CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
                 const std::string             a_realm,
                 const phase::which_phase      a_phase,
                 const int                     m_ncomp);
      ~CdrStorage();

      void
      allocateStorage(const int a_p);
      void
      deallocateStorage();

      EBAMRCellData&
      getScratch()
      {
        return m_scratch;
      }
      EBAMRCellData&
      getScratch2()
      {
        return m_scratch2;
      }
      EBAMRCellData&
      getError()
      {
        return m_error;
      }
      EBAMRCellData&
      getGradient()
      {
        return m_gradient;
      }
      EBAMRCellData&
      getOld()
      {
        return m_old;
      }
      EBAMRCellData&
      getScratchD()
      {
        return m_scratchD;
      }
      EBAMRCellData&
      getDivF()
      {
        return m_divF;
      }

      Vector<EBAMRCellData>&
      getPhi()
      {
        return m_phi;
      }
      Vector<EBAMRCellData>&
      getFAR()
      {
        return m_FAR;
      }
      Vector<EBAMRCellData>&
      getFD()
      {
        return m_FD;
      }
      Vector<EBAMRCellData>&
      getF()
      {
        return m_F;
      }

      EBAMRIVData&
      getEbScratch1()
      {
        return m_scratchIVs;
      }
      EBAMRIVData&
      getEbScratchD()
      {
        return m_scratchIVD;
      }
      EBAMRIVData&
      getEbState()
      {
        return m_scratchIV1;
      }
      EBAMRIVData&
      getEbVelo()
      {
        return m_scratchIV2;
      }
      EBAMRIVData&
      getEbFlux()
      {
        return m_scratchIV3;
      }
      EBAMRIVData&
      getEbGrad()
      {
        return m_scratchIV4;
      }

      EBAMRIFData&
      getDomainState()
      {
        return m_scratchIF1;
      }
      EBAMRIFData&
      getDomainVelo()
      {
        return m_scratchIF2;
      }
      EBAMRIFData&
      getDomainFlux()
      {
        return m_scratchIF3;
      }
      EBAMRIFData&
      getDomainGrad()
      {
        return m_scratchIF4;
      }

    protected:
      std::string            m_realm;
      RefCountedPtr<AmrMesh> m_amr;      // Needed for allocation. Will allocate
      phase::which_phase     m_phase;    // on this phase
      int                    m_ncomp;    // with this many components
      int                    m_p;        // IMEX_SDC order
      bool                   m_subcycle; // Extra storage for subcycling

      EBAMRCellData m_scratch;  // Scratch data
      EBAMRCellData m_scratch2; // Scratch data
      EBAMRCellData m_error;    // Error estimator
      EBAMRCellData m_gradient; // Storage for holding gradient
      EBAMRCellData m_old;      // "Old solution". Only used for subcycling
      EBAMRCellData m_scratchD; // SpaceDim scratch storage
      EBAMRCellData m_divF;     // SpaceDim scratch storage

      Vector<EBAMRCellData> m_phi; // phi^k in the  IMEX_SDC notation
      Vector<EBAMRCellData> m_FAR; // FAR(phi^k) in the IMEX_SDC notation
      Vector<EBAMRCellData> m_FD;  // FD(phi^k) in the IMEX_SDC notation
      Vector<EBAMRCellData> m_F;   // F(phi^k) in the IMEX_SDC notation

      EBAMRIVData m_scratchIVs; // Scratch data
      EBAMRIVData m_scratchIVD; // Scratch data
      EBAMRIVData m_scratchIV1; // Scratch data
      EBAMRIVData m_scratchIV2; // Scratch data
      EBAMRIVData m_scratchIV3; // Scratch data
      EBAMRIVData m_scratchIV4; // Scratch data

      EBAMRIFData m_scratchIF1; // Scratch data
      EBAMRIFData m_scratchIF2; // Scratch data
      EBAMRIFData m_scratchIF3; // Scratch data
      EBAMRIFData m_scratchIF4; // Scratch data
    };

    class CdrPlasmaImExSdcStepper::FieldStorage
    {
    public:
      FieldStorage();
      FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
                   const std::string             a_realm,
                   const phase::which_phase      a_phase,
                   const int                     a_ncomp);
      ~FieldStorage();

      void
      allocateStorage(const int a_p);
      void
      deallocateStorage();

      MFAMRCellData&
      getPrevious()
      {
        return m_previous;
      }
      EBAMRCellData&
      getElectricFieldCell()
      {
        return m_E_cell;
      }
      EBAMRFluxData&
      getElectricFieldFace()
      {
        return m_E_face;
      }
      EBAMRIVData&
      getElectricFieldEb()
      {
        return m_E_eb;
      }
      EBAMRIFData&
      getElectricFieldDomain()
      {
        return m_E_dom;
      }

    protected:
      std::string            m_realm;
      RefCountedPtr<AmrMesh> m_amr;   // Needed for allocation. Will allocate
      phase::which_phase     m_phase; // on this phase
      int                    m_ncomp; // with this many components
      int                    m_p;     // IMEX_SDC order

      MFAMRCellData m_previous; // Solution at previous time step
      EBAMRCellData m_E_cell;   // Cell-centered E
      EBAMRFluxData m_E_face;   // Face-centered E
      EBAMRIVData   m_E_eb;     // EB-centered E
      EBAMRIFData   m_E_dom;    // E-field on domain faces
    };

    class CdrPlasmaImExSdcStepper::RtStorage
    {
    public:
      RtStorage();
      RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
                const std::string             a_realm,
                const phase::which_phase      a_phase,
                const int                     a_ncomp);
      ~RtStorage();

      void
      allocateStorage(const int a_p);
      void
      deallocateStorage();

      EBAMRCellData&
      getPrevious()
      {
        return m_previous;
      }
      EBAMRIVData&
      getEbFlux()
      {
        return m_scratchIV;
      }
      EBAMRIFData&
      getDomainFlux()
      {
        return m_scratchIF;
      }

    protected:
      std::string            m_realm;
      RefCountedPtr<AmrMesh> m_amr;   // Needed for allocation. Will allocate
      phase::which_phase     m_phase; // on this phase
      int                    m_ncomp; // with this many components
      int                    m_p;     // IMEX_SDC order

      EBAMRCellData m_previous;  // Solution at previous time step
      EBAMRIVData   m_scratchIV; // Used for EB-centered isotropic density
      EBAMRIFData   m_scratchIF; // Used for domain-centered isotropic density
    };

    class CdrPlasmaImExSdcStepper::SigmaStorage
    {
    public:
      SigmaStorage();
      SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
                   const std::string             a_realm,
                   const phase::which_phase      a_phase,
                   const int                     a_ncomp);
      ~SigmaStorage();

      void
      allocateStorage(const int a_p);
      void
      deallocateStorage();

      EBAMRIVData&
      getScratch()
      {
        return m_scratch;
      }
      EBAMRIVData&
      getError()
      {
        return m_error;
      }

      Vector<EBAMRIVData>&
      getSigmaSolver()
      {
        return m_sigma;
      }
      Vector<EBAMRIVData>&
      getFnew()
      {
        return m_Fnew;
      }
      Vector<EBAMRIVData>&
      getFold()
      {
        return m_Fold;
      }

    protected:
      std::string            m_realm;
      RefCountedPtr<AmrMesh> m_amr;   // Needed for allocation. Will allocate
      phase::which_phase     m_phase; // on this phase
      int                    m_ncomp; // with this many components
      int                    m_p;     // IMEX_SDC order

      EBAMRIVData m_scratch; // Scratch data
      EBAMRIVData m_error;   // Error

      Vector<EBAMRIVData> m_sigma; // sigma^k in the IMEX_SDC notation
      Vector<EBAMRIVData> m_Fnew;  // New operator lsopes
      Vector<EBAMRIVData> m_Fold;  // Lagged operator slopes
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
