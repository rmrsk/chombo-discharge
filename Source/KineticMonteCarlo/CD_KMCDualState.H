/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCDualState.H
  @brief  Declaration of a simple state vector for running Kinetic Monte Carlo for plasma problems. 
  @author Robert Marskar
*/

#ifndef CD_KMCDualState_H
#define CD_KMCDualState_H

// Std includes
#include <vector>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a "dual state" for advancing with the Kinetic Monte Carlo module. 
  @details This state consists of both reactive and non-reactive species. The reactive species can appear on the left- and 
  right-hand side of the reaction. The non-reactive species can only appear on the right-hand side of the reaction. 
  @note The template parameter indicates the integer type used for tracking the states. 
*/
template <typename T = long long>
class KMCDualState
{
public:
  using State = std::vector<T>;

  /*!
    @brief Default constructor. 
  */
  inline KMCDualState() = default;

  /*!
    @brief Copy constructor
  */
  inline KMCDualState(const KMCDualState&) = default;

  /*!
    @brief Disallowed move constructor
  */
  inline KMCDualState(const KMCDualState&&) = delete;

  /*!
    @brief Full constructor. 
    @param[in] a_numReactiveSpecies    Number of reactive species in the state vector.
    @param[in] a_numNonReactiveSpecies Number of non-reactive species in the state vector.
  */
  inline KMCDualState(const size_t a_numReactiveSpecies, const size_t a_numNonReactiveSpecies) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~KMCDualState();

  /*!
    @brief Copy assignment constructor
  */
  inline KMCDualState&
  operator=(const KMCDualState&) = default;

  /*!
    @brief Disallowed move assignement
  */
  inline KMCDualState&
  operator=(const KMCDualState&&) = delete;

  /*!
    @brief Define function constructor. 
    @param[in] a_numReactiveSpecies    Number of reactive species in the state vector.
    @param[in] a_numNonReactiveSpecies Number of non-reactive species in the state vector.
  */
  inline void
  define(const size_t a_numReactiveSpecies, const size_t a_numNonReactiveSpecies) noexcept;

  /*!
    @brief Check if state is a valid state. An invalid state will have a negative number of reactants/non-reactants. 
  */
  inline bool
  isValidState() const noexcept;

  /*!
    @brief Get modifiable reactive state. 
    @return m_reactiveState
  */
  inline State&
  getReactiveState() noexcept;

  /*!
    @brief Get reactive state. 
    @return m_reactiveState
  */
  inline const State&
  getReactiveState() const noexcept;

  /*!
    @brief Get modifiable non-reactive state
    @return m_nonReactiveState
  */
  inline State&
  getNonReactiveState() noexcept;

  /*!
    @brief Get non-reactive state
    @return m_nonReactiveState
  */
  inline const State&
  getNonReactiveState() const noexcept;

  /*!
    @brief Linearize the state onto an output vector.
    @details This is mostly by KMCSolver when computing the necessary Jacobians for implicit updates
  */
  inline std::vector<T>
  linearOut() const noexcept;

  /*!
    @brief Linearize the input buffer.
    @details This routine should take the input vector and linearize it onto the current state. This is mostly used
    by KMCSolver when computing the Jacobians for implicit updates.
    @param[in] a_linearizedState Linearized state vector
  */
  inline void
  linearIn(const std::vector<T>& a_linearizedState) noexcept;

protected:
  /*!
    @brief Number of reactive species
  */
  int m_numReactiveSpecies;

  /*!
    @brief Number of non-reactive species
  */
  int m_numNonReactiveSpecies;

  /*!
    @brief Reactive state.
  */
  State m_reactiveState;

  /*!
    @brief Non-reactive state. 
  */
  State m_nonReactiveState;
};

/*!
  @brief State printing function
  @param[in] ostr    Output stream
  @param[in] a_state State vector
*/
template <typename T>
inline std::ostream&
operator<<(std::ostream& ostr, const KMCDualState<T>& a_state);

#include <CD_NamespaceFooter.H>

#include <CD_KMCDualStateImplem.H>

#endif
