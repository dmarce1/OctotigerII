/** @file
 * @brief Compile-time CGS dimensions and quantity-safe arithmetic.
 * @ingroup units
 */
#pragma once
#include <algorithm>
#include <boost/units/base_units/cgs/centimeter.hpp>
#include <boost/units/base_units/cgs/gram.hpp>
#include <boost/units/base_units/si/kelvin.hpp>
#include <boost/units/base_units/si/second.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/derived_dimension.hpp>
#include <boost/units/make_system.hpp>
#include <boost/units/quantity.hpp>
#include <cmath>
#include "octotigerII/math/Real.hpp"


namespace octotigerII::units {
// CGS mechanical units with kelvin for thermodynamics. No SI mechanical
// quantities or selectable code-unit scales are used by the application.
using System = boost::units::make_system<boost::units::cgs::centimeter_base_unit, boost::units::cgs::gram_base_unit, boost::units::si::second_base_unit,
	boost::units::si::kelvin_base_unit>::type;
template <typename Dimension, int Exponent, typename Tail>
using DimensionTerm = std::conditional_t<Exponent == 0, Tail, boost::units::list<boost::units::dim<Dimension, boost::units::static_rational<Exponent>>, Tail>>;
template <int L, int M, int T, int K = 0>
using Quantity = boost::units::quantity<
	boost::units::unit<typename boost::units::make_dimension_list<DimensionTerm<boost::units::length_base_dimension, L,
						   DimensionTerm<boost::units::mass_base_dimension, M,
							   DimensionTerm<boost::units::time_base_dimension, T,
								   DimensionTerm<boost::units::temperature_base_dimension, K, boost::units::dimensionless_type>>>>>::type,
		System>,
	Real>;
using Length = Quantity<1, 0, 0>;
using Mass = Quantity<0, 1, 0>;
using Time = Quantity<0, 0, 1>;
using Temperature = Quantity<0, 0, 0, 1>;
using Velocity = Quantity<1, 0, -1>;
using VelocitySquared = Quantity<2, 0, -2>;
using Acceleration = Quantity<1, 0, -2>;
using Density = Quantity<-3, 1, 0>;
using MomentumDensity = Quantity<-2, 1, -1>;
using Pressure = Quantity<-1, 1, -2>;
using EnergyDensity = Pressure;
using Energy = Quantity<2, 1, -2>;
using EnergyFlux = Quantity<0, 1, -3>;
using EnergyFluxTransport = Quantity<1, 1, -4>;
using MassFlux = MomentumDensity;
using MomentumFlux = Pressure;
using Volume = Quantity<3, 0, 0>;
using Momentum = Quantity<1, 1, -1>;
using InverseTime = Quantity<0, 0, -1>;
using MassPerLength = Quantity<-1, 1, 0>;
using TimePerLength = Quantity<-1, 0, 1>;
using GravitationalConstant = Quantity<3, -1, -2>;
using BoltzmannConstant = Quantity<2, 1, -2, -1>;
using RadiationConstant = Quantity<-1, 1, -2, -4>;
using Action = Quantity<2, 1, -1>;

template <typename U>
constexpr Real value(boost::units::quantity<U, Real> const& q) {
	return q.value();
}

template <typename U>
bool finite(boost::units::quantity<U, Real> const& q) {
	using std::isfinite;

	return isfinite(q.value());
}

template <typename Q>
Q hypot(Q a, Q b) {
	using std::sqrt;

	auto scale = std::max(boost::units::abs(a), boost::units::abs(b));
	if (scale == Q{}) return Q{};
	Real x = a / scale, y = b / scale;
	return scale * sqrt(x * x + y * y);
}
using boost::units::abs;
using boost::units::sqrt;


}	 // namespace octotigerII::units


// ADL serialization for HPX archives. Quantities are reconstructed with the
// same compile-time unit; the wire format contains only their CGS value.
namespace boost::units {
/// Serialize this value with its compile-time quantity types preserved.
template <typename Archive, typename Unit>
void serialize(Archive& archive, quantity<Unit, double>& q, unsigned) {
	double v = q.value();
	archive & v;
	q = quantity<Unit, double>::from_value(v);
}
}	 // namespace boost::units
