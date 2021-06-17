/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GeoCoarsener.H
  @brief  Declaration of a utility class for removing geometrically tagged cells
  @author Robert Marskar
*/

#ifndef CD_GeoCoarsener_H
#define CD_GeoCoarsener_H

// Chombo includes
#include <Vector.H>
#include <IntVectSet.H>

// Our includes
#include <CD_RealBox.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Utility class for removing geometrically tagged cells
  @details This class is an optimization class for coarsening of EBs. The usage of this class is simple; the user defines m_coarsen_boxes and a level for each coarsening box. These are passed to plasma_engine which will remove geometric tags in those regions (above the level specified)
*/
class GeoCoarsener {
public:

  /*!
    @brief Constructor
  */
  GeoCoarsener();

  /*!
    @brief Destructor
  */
  virtual ~GeoCoarsener();

  /*!
    @brief Coarsen tags
    @param[in] a_tags   Tags to be modified. 
    @param[in] a_dx     Grid resolutions. 
    @param[in] a_probLo Domain origin. 
  */
  virtual void coarsenTags(Vector<IntVectSet>& a_tags, const Vector<Real>& a_dx, const RealVect& a_probLo) const;

  /*!
    @brief Get coarsening boxes
  */
  virtual Vector<RealBox> getCoarsenBoxes();

  /*!
    @brief Get coarsening levels
  */
  virtual Vector<int> getCoarsenLevels();
  
protected:

  /*!
    @brief Coarsening boxes
  */
  Vector<RealBox> m_coarsen_boxes;

  /*!
    @brief Coarsening levels. Tags are removed above this level
  */
  Vector<int> m_coarsen_levels;

  /*!
    @brief Inverse box or not
  */
  Vector<int> m_inverse;
};

#include <CD_NamespaceFooter.H>

#endif