/*!
  @file   point_mass.cpp
  @brief  Implementation of point_mass.H
  @author Robert Marskar
  @date   May 2020
*/

#include "point_mass.H"

point_mass::point_mass(){

}

point_mass::point_mass(const RealVect a_pos, const Real a_mass){
  m_pos  = a_pos;
  m_mass = a_mass;
}

point_mass::point_mass(const std::vector<point_mass>& a_point_masses){
  m_pos  = RealVect::Zero;
  m_mass = 0.0;

  for (int i = 0; i < a_point_masses.size(); i++){
    const RealVect& p = a_point_masses[i].m_pos;
    const Real& m     = a_point_masses[i].m_mass;
    
    m_mass += m;
    m_pos  += m*p;
  }

  m_pos = m_pos/m_mass;
}

point_mass::~point_mass(){

}