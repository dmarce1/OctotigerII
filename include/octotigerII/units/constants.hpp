/** @file
 * @brief Physical constants in CGS from the 2022 CODATA adjustment.
 * @ingroup units
 */
#pragma once
#include "octotigerII/units/cgs.hpp"


namespace octotigerII::constants {
// CODATA 2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
/// Constants use the 2022 adjustment: @ref ref_mohr2025 "Mohr et al. (2025)".
/// The numerical literal for G is in cm³/(g s²).
inline constexpr auto G = units::GravitationalConstant::from_value(6.67430e-8);			 // cm^3/(g s^2)
inline constexpr auto c = units::Velocity::from_value(2.99792458e10);					 // cm/s, exact
inline constexpr auto atomicMassUnit = units::Mass::from_value(1.66053906892e-24);		 // g
inline constexpr auto boltzmann = units::BoltzmannConstant::from_value(1.380649e-16);	 // erg/K, exact
inline constexpr auto planck = units::Action::from_value(6.62607015e-27);				 // erg s, exact
// a = 8*pi^5*k_B^4/(15*h^3*c^3), rounded to double precision.
inline constexpr auto radiation = units::RadiationConstant::from_value(7.565733250280004e-15);	  // erg/(cm^3 K^4)
inline constexpr auto m_u = atomicMassUnit;
inline constexpr auto k_B = boltzmann;
inline constexpr auto a_r = radiation;
}	 // namespace octotigerII::constants
