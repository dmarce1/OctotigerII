/** @file
 * @brief Conservative, local gray gas--radiation exchange in inertial variables.
 */
#pragma once
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/radiation/radiationTransport.hpp"

namespace octotigerII::radiation {

/// Gray absorption opacity per unit mass, in cm^2/g. Scattering is not included.
using Opacity = units::Quantity<2, -1, 0>;

/// Integrate the equal-opacity first-order mixed-frame terms from
/// Skinner--Ostriker (2013), equations (5), retaining the full M1 tensor.
/// With sigma=rho*opacity, q=sigma[c(E-aT^4)-v.F/c] and
/// G=sigma/c[F-P.v-aT^4 v], the gas receives (q,G) and radiation receives
/// (-ratio*q,-c^2*ratio*G). Velocities and moments are inertial, including on
/// a rotating mesh. The conserved source combinations are Egas+E/ratio and
/// momentum+F/(c^2*ratio); only ratio=1 gives the physical total invariants.
/// The additional second-order terms in equations (5) are not included.
///
/// A two-stage L-stable SDIRK method is second order for resolved smooth source
/// evolution. Newton stages and accepted steps must have positive gas thermal
/// energy and realizable radiation to scaled floating-point roundoff. Failed
/// steps are subdivided, never clipped.
/// The selected dual-energy thermal state is advanced by gas total-energy change
/// minus kinetic-energy change, and its auxiliary is explicitly updated.
/// Both arguments are unchanged on failure. Zero opacity/interval is a no-op.
void couple(hydro::ConservedState& gas, RadiationSystem::State& radiation,
	hydro::HydroSystem const& system, Opacity opacity, Real lightSpeedRatio, units::Time interval);

/// Integrate the same sources together with constant transport/gravity driving.
/// The increments are the accepted finite-volume changes over interval, not
/// rates. They are included exactly once in the returned states. Each implicit
/// stage includes the radiation pressure-divergence drive, preserving a stiff
/// diffusion balance that a transport/source splitting would destroy.
/// Gas density follows its supplied increment. The selected-thermal discrepancy
/// from Egas-kinetic is interpolated between the uncoupled gas endpoints.
void coupleForced(hydro::ConservedState& gas, RadiationSystem::State& radiation,
	hydro::ConservedState const& gasIncrement, RadiationSystem::State const& radiationIncrement,
	hydro::HydroSystem const& system, Opacity opacity, Real lightSpeedRatio, units::Time interval);

} // namespace octotigerII::radiation
