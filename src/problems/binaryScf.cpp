#include "octotigerII/problems/binaryScf.hpp"
#include "octotigerII/problems/bipolytrope.hpp"
#include "octotigerII/gravity/isolatedPotential.hpp"
#include "octotigerII/profiling.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <limits>
#include <mutex>
#include <numeric>
#include <numbers>
#include <sstream>
#include <stdexcept>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/synchronization/mutex.hpp>
#endif

namespace octotigerII::problems {
namespace {
	using Pair = std::array<Real, 2>;
	using Vector = std::array<Real, 3>;
	std::size_t index(int x, int y, int z, int n) { return (std::size_t(z) * n + y) * n + x; }
	// The synchronous binary is reflection symmetric about its orbital plane
	// and the plane through both centers and the rotation axis. Project out
	// transverse translation modes; otherwise roundoff plus acceleration can
	// move a star away from the axis assumed by the Roche-saddle search.
	void reflectSymmetry(std::vector<Real>& a, int n) {
		for (int z = 0; z < n/2; ++z) for (int y = 0; y < n/2; ++y) for (int x = 0; x < n; ++x) {
			std::array<std::size_t, 4> const ids{index(x,y,z,n),index(x,n-1-y,z,n),index(x,y,n-1-z,n),index(x,n-1-y,n-1-z,n)};
			Real const value = .25*((a[ids[0]]+a[ids[1]])+(a[ids[2]]+a[ids[3]]));
			for (auto i : ids) a[i] = value;
		}
	}
	// Keep only the connected stellar lobe. The centrifugal potential also
	// admits an unphysical exterior branch; clipping H alone would populate it.
	void retainLobe(std::vector<Real>& rho, std::size_t peak, int n) {
		std::vector<unsigned char> visited(rho.size(), 0);
		std::vector<std::size_t> queue{peak};
		visited[peak] = 1;
		for (std::size_t head = 0; head < queue.size(); ++head) {
			auto const i = queue[head];
			std::array<int, 3> const c{int(i % n), int(i / n % n), int(i / (n * n))};
			for (int d = 0; d < 3; ++d) {
				if (c[d] == 0 || c[d] == n - 1)
					throw std::runtime_error("SCF stellar lobe reaches reference boundary; requested binary is unbounded or unresolved");
				std::ptrdiff_t const stride = d == 0 ? 1 : d == 1 ? n : n * n;
				for (int sign : {-1, 1}) {
					auto const j = std::size_t(std::ptrdiff_t(i) + sign * stride);
					if (!visited[j] && rho[j] > 0) { visited[j] = 1; queue.push_back(j); }
				}
			}
		}
		for (std::size_t i = 0; i < rho.size(); ++i) if (!visited[i]) rho[i] = 0;
	}

	/// Bounded-history Anderson mixing of the normalized fixed-point map.
	/// The acceptance test still uses the original, undamped density correction.
	class Mixing {
		std::vector<std::vector<Real>> errors_, outputs_;
		Real previous_ = std::numeric_limits<Real>::infinity();
	public:
		std::vector<Real> operator()(std::vector<Real> const& x, std::vector<Real> const& f,
			Real relaxation, int history, Real residual) {
			if (residual > 2*previous_) { errors_.clear(); outputs_.clear(); }
			previous_ = residual;
			std::vector<Real> error(x.size()), output(x.size());
			for (std::size_t i = 0; i < x.size(); ++i) { error[i] = f[i]-x[i]; output[i] = x[i]+relaxation*error[i]; }
			errors_.push_back(std::move(error)); outputs_.push_back(output);
			if (int(errors_.size()) > history+1) { errors_.erase(errors_.begin()); outputs_.erase(outputs_.begin()); }
			int const m = int(errors_.size())-1;
			std::vector<std::vector<Real>> delta(m, std::vector<Real>(x.size()));
			for (int j = 0; j < m; ++j) for (std::size_t i = 0; i < x.size(); ++i) delta[j][i] = errors_[j+1][i]-errors_[j][i];
			std::vector<std::vector<Real>> matrix(m, std::vector<Real>(m+1));
			Real scale = 0;
			for (int j = 0; j < m; ++j) {
				for (int k = 0; k < m; ++k) matrix[j][k] = std::inner_product(delta[j].begin(),delta[j].end(),delta[k].begin(),Real(0));
				matrix[j][m] = std::inner_product(delta[j].begin(),delta[j].end(),errors_.back().begin(),Real(0));
				scale = std::max(scale,matrix[j][j]);
			}
			if (!(scale > 0)) return output;
			for (int j = 0; j < m; ++j) matrix[j][j] += 1e-10*scale;
			for (int j = 0; j < m; ++j) {
				int pivot = j;
				for (int k = j+1; k < m; ++k) if (std::abs(matrix[k][j]) > std::abs(matrix[pivot][j])) pivot = k;
				std::swap(matrix[j],matrix[pivot]);
				Real const divisor = matrix[j][j];
				if (std::abs(divisor) < 1e-14*scale) return output;
				for (int k = j; k <= m; ++k) matrix[j][k] /= divisor;
				for (int row = 0; row < m; ++row) if (row != j) {
					Real const factor = matrix[row][j];
					for (int k = j; k <= m; ++k) matrix[row][k] -= factor*matrix[j][k];
				}
			}
			for (int j = 0; j < m; ++j) {
				if (!std::isfinite(matrix[j][m]) || std::abs(matrix[j][m]) > 100) return outputs_.back();
				for (std::size_t i = 0; i < x.size(); ++i) output[i] -= matrix[j][m]*(outputs_[j+1][i]-outputs_[j][i]);
			}
			// A large extrapolation or appreciable negative mass uses simple mixing.
			Real change = 0, negative = 0, mass = 0;
			for (std::size_t i = 0; i < x.size(); ++i) { change += std::abs(output[i]-x[i]); negative += std::max(Real(0),-output[i]); mass += x[i]; }
			if (change > .5*mass || negative > .001*mass) return outputs_.back();
			for (auto& r : output) r = std::max(Real(0),r);
			return output;
		}
	};
}

void BinaryScf::validate(Config const& c) {
	static_assert(ndim == 3);
	auto const& p = c.scf;
	if (c.hydro.eos != "ideal" || std::abs(c.hydro.gamma - Real(5) / 3) > 1e-12 || c.radiationEnabled() || c.hasExternalAcceleration())
		throw std::invalid_argument("binary-scf requires ideal-gas gamma=5/3, gravity, and no radiation or external acceleration");
	if (!c.mesh.boundary.all(finiteVolume::BoundaryCondition::Outflow))
		throw std::invalid_argument("binary-scf requires isolated outflow boundaries");
	if (!(p.primaryMass > units::Mass{}) || !units::finite(p.primaryMass) ||
		!(p.separation > units::Length{}) || !units::finite(p.separation))
		throw std::invalid_argument("SCF masses and separation must be finite and positive");
	for (auto v : {p.massRatio, p.atmosphereFraction, p.referenceWidth, p.tolerance, p.relaxation, p.virialTolerance})
		if (!(v > 0) || !std::isfinite(v)) throw std::invalid_argument("SCF controls must be finite and positive");
	if (p.referenceWidth < 1.5 || p.referenceWidth > 10)
		throw std::invalid_argument("SCF referenceWidth must lie in [1.5,10]");
	if (!(p.evolveOrbits >= 0) || !std::isfinite(p.evolveOrbits) || p.framesPerOrbit < 0 || p.framesPerOrbit > 10000)
		throw std::invalid_argument("SCF requires nonnegative finite evolveOrbits and framesPerOrbit=0..10000");
	if (p.cells < 16 || p.cells > 128 || (p.cells & (p.cells - 1)) || p.maxIterations < 1 || p.history < 0 || p.history > 6 ||
		p.tolerance >= 0.1 || p.relaxation > 1 || p.atmosphereFraction >= 1e-3)
		throw std::invalid_argument("SCF requires power-of-two cells=16..128, iterations>0, tolerance<0.1, relaxation<=1, atmosphereFraction<1e-3");
	for (int s = 0; s < 2; ++s) {
		Bipolytrope const eos(p.coreIndex[s], p.envelopeIndex[s], p.interfaceFraction[s], p.densityJump[s], 1, 1);
		if (!(p.fill[s] > 0 && p.fill[s] <= 1) || !std::isfinite(p.fill[s]))
			throw std::invalid_argument("SCF fill must be in (0,1]; donor.fill=1 selects Roche-lobe contact");
	}
}

BinaryScf::BinaryScf(Config const& c) : parameters_(c.scf), n_(c.scf.cells), width_(0), lower_(-.5*c.scf.referenceWidth), atmosphere_(c.scf.atmosphereFraction) {
	profiling::Region scfProfile("scf.solve");
	validate(c);
	width_ = c.scf.referenceWidth/n_;
	auto const& p = c.scf;
	std::size_t const count = std::size_t(n_) * n_ * n_;
	Real const dv = width_ * width_ * width_;
	Pair const target{1 / (1 + p.massRatio), p.massRatio / (1 + p.massRatio)};
	std::array<std::vector<Real>, 2> rho{std::vector<Real>(count), std::vector<Real>(count)};
	std::vector<Vector> positions(count);
	for (int z = 0; z < n_; ++z) for (int y = 0; y < n_; ++y) for (int x = 0; x < n_; ++x) {
		auto const i = index(x, y, z, n_);
		positions[i] = {lower_ + (x + .5) * width_, lower_ + (y + .5) * width_, lower_ + (z + .5) * width_};
		for (int s = 0; s < 2; ++s) {
			Real const q = target[s] / target[1-s], q13 = std::cbrt(q);
			Real const radius = 0.49 * q13 * q13 / (0.6 * q13 * q13 + std::log1p(q13)) * p.fill[s];
			Real const center = s == 0 ? -target[1] : target[0];
			Real const r = std::hypot(positions[i][0] - center, positions[i][1], positions[i][2]) / radius;
			rho[s][i] = r < 1 ? std::pow(1 - r*r, 1.5) : 0;
		}
	}
	Real orbitalJ = 0;
	std::vector<Real> phi(count), effective(count);
	std::array<std::unique_ptr<Bipolytrope>, 2> eos;
	std::array<std::size_t, 2> peaks{};
	Real residual = 1;
	bool converged = false;
	gravity::IsolatedPotential const potential(n_, width_);
	Mixing mixing;
	for (int iteration = 0; iteration <= p.maxIterations; ++iteration) {
		// Enforce target masses before solving gravity; the undamped residual below
		// also tests this normalization, so it cannot hide a persistent mass error.
		Pair mass{}, center{}, peakDensity{};
		for (int s = 0; s < 2; ++s) {
			reflectSymmetry(rho[s],n_);
			mass[s] = std::accumulate(rho[s].begin(), rho[s].end(), Real(0)) * dv;
			if (!(mass[s] > 0) || !std::isfinite(mass[s])) throw std::runtime_error("SCF lost a star; increase resolution or change parameters");
			for (std::size_t i = 0; i < count; ++i) {
				rho[s][i] *= target[s] / mass[s];
				center[s] += positions[i][0] * rho[s][i] * dv / target[s];
				if (rho[s][i] > peakDensity[s]) { peakDensity[s] = rho[s][i]; peaks[s] = i; }
			}
		}
		Real const axis = target[0] * center[0] + target[1] * center[1];
		Real const separation = center[1] - center[0];
		if (!(separation > 4 * width_)) throw std::runtime_error("SCF binary centers merged or are unresolved");
		Real const inertia = target[0] * target[1] * separation * separation;
		if (iteration == 0) orbitalJ = inertia / std::pow(separation, 1.5);
		Real const omega = orbitalJ / inertia;
		std::vector<Real> sources(count);
		for (std::size_t i = 0; i < count; ++i) sources[i] = rho[0][i] + rho[1][i];
		{
			profiling::Region gravityProfile("scf.gravity_fft");
			phi = potential(sources);
		}
		reflectSymmetry(phi,n_);
		for (std::size_t i = 0; i < count; ++i) {
			effective[i] = phi[i] - .5 * omega * omega * (std::pow(positions[i][0] - axis, 2) + positions[i][1] * positions[i][1]);
		}
		// Independent point-mass evaluation on the symmetry axis avoids moving
		// the Roche saddle when the reference grid is translated by a fraction
		// of a cell. The central grid plane contains no source-cell centers.
		Real l1, l1phi;
		{
		profiling::Region saddleProfile("scf.l1_search");
		auto axisPotential = [&](Real x) {
			Real value = -.5*omega*omega*(x-axis)*(x-axis);
			for (std::size_t i = 0; i < count; ++i) if (sources[i] > 0)
				value -= sources[i]*dv/std::hypot(x-positions[i][0],positions[i][1],positions[i][2]);
			return value;
		};
		Real lo = center[0], hi = center[1];
		constexpr Real ratio = .6180339887498948482;
		Real left = hi-ratio*(hi-lo), right = lo+ratio*(hi-lo);
		Real fl = axisPotential(left), fr = axisPotential(right);
		for (int j = 0; j < 60; ++j) {
			if (fl < fr) { lo = left; left = right; fl = fr; right = lo+ratio*(hi-lo); fr = axisPotential(right); }
			else { hi = right; right = left; fr = fl; left = hi-ratio*(hi-lo); fl = axisPotential(left); }
		}
		l1 = (lo+hi)/2; l1phi = axisPotential(l1);
		}
		Pair minimum{std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity()}, bernoulli{};
		for (std::size_t i = 0; i < count; ++i) {
			int const s = positions[i][0] < l1 ? 0 : 1;
			// Restrict the minimum search to the current star, excluding the
			// unbounded centrifugal exterior.
			if (rho[s][i] > peakDensity[s] * .01 && effective[i] < minimum[s]) { minimum[s] = effective[i]; peaks[s] = i; }
		}
		for (int s = 0; s < 2; ++s) {
			bernoulli[s] = minimum[s] + p.fill[s] * (l1phi - minimum[s]);
			if (!(bernoulli[s] > minimum[s])) throw std::runtime_error("SCF has no resolved inner potential saddle");
			eos[s] = std::make_unique<Bipolytrope>(p.coreIndex[s], p.envelopeIndex[s], p.interfaceFraction[s], p.densityJump[s],
				peakDensity[s], bernoulli[s] - minimum[s]);
		}
		auto next = rho;
		{
			profiling::Region densityProfile("scf.density_update");
			for (std::size_t i = 0; i < count; ++i) for (int s = 0; s < 2; ++s)
				next[s][i] = ((positions[i][0] < l1) == (s == 0)) ? eos[s]->updateDensity(bernoulli[s] - effective[i],rho[s][i]) : 0;
			for (int s = 0; s < 2; ++s) retainLobe(next[s], peaks[s], n_);
		}
		residual = 0;
		Real enthalpyResidual = 0;
		for (int s = 0; s < 2; ++s) {
			Real error = 0, enthalpyError = 0, enthalpyNorm = 0;
			for (std::size_t i = 0; i < count; ++i) {
				error += std::abs(next[s][i] - rho[s][i]) * dv;
				enthalpyError += rho[s][i]*std::abs(eos[s]->enthalpy(rho[s][i])+effective[i]-bernoulli[s])*dv;
				enthalpyNorm += rho[s][i]*(bernoulli[s]-minimum[s])*dv;
			}
			if (!std::isfinite(error) || !std::isfinite(enthalpyError) || !(enthalpyNorm > 0) || !std::isfinite(enthalpyNorm))
				throw std::runtime_error("SCF generated a nonfinite or degenerate stellar residual");
			residual = std::max(residual, error / target[s]);
			enthalpyResidual = std::max(enthalpyResidual,enthalpyError/enthalpyNorm);
		}
		diagnostics_.iterations = iteration;
		diagnostics_.densityResidual = residual;
		diagnostics_.maxStarBernoulliResidual = enthalpyResidual;
		diagnostics_.mass = target;
		diagnostics_.center = center;
		diagnostics_.centralDensity = peakDensity;
		diagnostics_.rotationCenter = axis;
		diagnostics_.separation = separation;
		diagnostics_.omega = omega;
		diagnostics_.l1 = l1;
		diagnostics_.l1Potential = l1phi;
		diagnostics_.bernoulli = bernoulli;
		diagnostics_.orbitalAngularMomentum = orbitalJ;
		if (iteration % 20 == 0) std::clog << "SCF iteration=" << iteration << " densityResidual=" << residual << " enthalpyResidual=" << enthalpyResidual << " separation=" << separation << " omega=" << omega << '\n';
		if (!std::isfinite(residual) || !std::isfinite(enthalpyResidual)) throw std::runtime_error("SCF generated a nonfinite residual");
		if (residual < p.tolerance && enthalpyResidual < p.tolerance) { converged = true; break; }
		if (iteration == p.maxIterations) break;
		profiling::Region mixingProfile("scf.mixing");
		std::vector<Real> x(2*count), f(2*count);
		for (int s = 0; s < 2; ++s) {
			Real const nextMass = std::accumulate(next[s].begin(), next[s].end(), Real(0))*dv;
			for (std::size_t i = 0; i < count; ++i) { x[s*count+i] = rho[s][i]/target[s]; f[s*count+i] = next[s][i]/nextMass; }
		}
		auto const update = mixing(x,f,p.relaxation,p.history,residual);
		for (int s = 0; s < 2; ++s) for (std::size_t i = 0; i < count; ++i) rho[s][i] = update[s*count+i]*target[s];
	}
	if (!converged) {
		std::ostringstream message;
		message << "SCF did not converge: density residual=" << residual << ", maximum stellar Bernoulli residual="
			<< diagnostics_.maxStarBernoulliResidual << "; increase scf.maxIterations/resolution or adjust relaxation/structure";
		throw std::runtime_error(message.str());
	}
	// The published density and its gravity belong to the same final iteration.
	// No density/energy change follows this diagnostic pass.
	cells_.resize(count);
	Real hError = 0, hNorm = 0;
	for (std::size_t i = 0; i < count; ++i) {
		auto& cell = cells_[i];
		for (int s = 0; s < 2; ++s) {
			Real const r = rho[s][i], core = r * eos[s]->coreFraction(r);
			if (r > 0) diagnostics_.volume[s] += dv;
			cell.density += r;
			cell.pressure += eos[s]->pressure(r);
			cell.partial[2*s] = core;
			cell.partial[2*s+1] = r - core;
			diagnostics_.coreMass[s] += core * dv;
			Real const h = eos[s]->enthalpy(r);
			hError += r * std::abs(h + effective[i] - diagnostics_.bernoulli[s]) * dv;
			hNorm += r * (diagnostics_.bernoulli[s] - effective[peaks[s]]) * dv;
			Real const spinRadius2 = std::pow(positions[i][0] - diagnostics_.center[s], 2) + positions[i][1] * positions[i][1];
			diagnostics_.spinAngularMomentum += r * spinRadius2 * diagnostics_.omega * dv;
		}
		Real const r2 = std::pow(positions[i][0] - diagnostics_.rotationCenter, 2) + positions[i][1] * positions[i][1];
		diagnostics_.kinetic += .5 * cell.density * diagnostics_.omega * diagnostics_.omega * r2 * dv;
		diagnostics_.potential += .5 * cell.density * phi[i] * dv;
		diagnostics_.pressureIntegral += cell.pressure * dv;
	}
	diagnostics_.bernoulliResidual = hError / hNorm;
	diagnostics_.virialResidual = std::abs(2 * diagnostics_.kinetic + diagnostics_.potential + 3 * diagnostics_.pressureIntegral) / std::abs(diagnostics_.potential);
	if (!std::isfinite(diagnostics_.virialResidual) || diagnostics_.virialResidual > p.virialTolerance)
		throw std::runtime_error("SCF density converged but virial residual=" + std::to_string(diagnostics_.virialResidual)
			+ " exceeds scf.virialTolerance; increase scf.cells or explicitly loosen the gate for a resolution study");
	length_ = p.separation / diagnostics_.separation;
	density_ = (p.primaryMass * (1 + p.massRatio)) / (length_ * length_ * length_);
	pressure_ = constants::G * density_ * density_ * length_ * length_;
	if (!units::finite(length_) || !units::finite(density_) || !units::finite(pressure_) || !units::finite(angularVelocity()) ||
		!(length_ > units::Length{}) || !(density_ > units::Density{}) || !(pressure_ > units::Pressure{}))
		throw std::invalid_argument("SCF physical scaling is outside the representable CGS range");
	if (atmosphere_ * density_ < units::Density::from_value(1e-14) || atmosphere_ * pressure_ < units::Pressure::from_value(1e-14))
		throw std::invalid_argument("SCF atmosphere lies below the hydro floors at the requested physical scale");
	std::clog << "SCF converged iterations=" << diagnostics_.iterations << " densityResidual=" << residual
		<< " bernoulliResidual=" << diagnostics_.bernoulliResidual << " virialResidual=" << diagnostics_.virialResidual
		<< " omega=" << units::value(angularVelocity()) << " rad/s\n";
}

std::shared_ptr<BinaryScf const> BinaryScf::get(Config const& c) {
	validate(c);
	auto const& p = c.scf;
	std::vector<Real> key{units::value(p.primaryMass), units::value(p.separation), p.massRatio, p.atmosphereFraction, p.referenceWidth,
		Real(p.cells), Real(p.maxIterations), Real(p.history), p.tolerance, p.relaxation, p.virialTolerance};
	for (auto const& pair : {p.coreIndex, p.envelopeIndex, p.interfaceFraction, p.densityJump, p.fill}) key.insert(key.end(), pair.begin(), pair.end());
#ifdef OCTOTIGERII_WITH_HPX
	static hpx::mutex mutex;
#else
	static std::mutex mutex;
#endif
	static std::map<std::vector<Real>, std::shared_ptr<BinaryScf const>> models;
	std::lock_guard lock(mutex);
	if (auto const found = models.find(key); found != models.end()) return found->second;
	auto result = std::make_shared<BinaryScf>(c);
	models.emplace(std::move(key), result);
	return result;
}

units::InverseTime BinaryScf::angularVelocity() const { return diagnostics_.omega * units::sqrt(constants::G * density_); }
units::Time BinaryScf::orbitalPeriod() const { return (2 * std::numbers::pi) / angularVelocity(); }

BinaryScf::Cell BinaryScf::sample(mesh::PhysicalCoordinates const& physical) const {
	Vector x{};
	for (int d = 0; d < 3; ++d) {
		if (!units::finite(physical[d])) throw std::invalid_argument("Nonfinite SCF sample coordinate");
		x[d] = Real(physical[d] / length_);
	}
	x[0] += diagnostics_.rotationCenter;
	std::array<int, 3> base{};
	Vector f{};
	for (int d = 0; d < 3; ++d) {
		Real const t = (x[d] - lower_) / width_ - .5;
		if (t < 0 || t > n_ - 1) return {};
		base[d] = std::min(n_ - 2, int(t)); f[d] = t - base[d];
	}
	Cell result;
	for (int k = 0; k < 8; ++k) {
		Real w = 1; auto c = base;
		for (int d = 0; d < 3; ++d) { bool const plus = (k >> d) & 1; w *= plus ? f[d] : 1 - f[d]; c[d] += plus; }
		auto const& a = cells_[index(c[0], c[1], c[2], n_)];
		result.density += w * a.density; result.pressure += w * a.pressure;
		for (int s = 0; s < 4; ++s) result.partial[s] += w * a.partial[s];
	}
	return result;
}

hydro::PrimitiveState BinaryScf::operator()(mesh::PhysicalCoordinates const& x) const {
	return state(sample(x), x);
}

hydro::PrimitiveState BinaryScf::state(Cell const& cell, mesh::PhysicalCoordinates const& x) const {
	hydro::PrimitiveState q;
	q.density() = (cell.density + atmosphere_) * density_;
	q.pressure() = (cell.pressure + atmosphere_) * pressure_;
	// A uniform, tenuous co-rotating atmosphere is explicit added mass/pressure.
	q.velocity(0) = -angularVelocity() * x[1];
	q.velocity(1) = angularVelocity() * x[0];
	return q;
}

BinaryScf::Cell BinaryScf::average(mesh::PhysicalCoordinates const& physical, units::Length cellWidth) const {
	Real const h = cellWidth / length_;
	if (!(h > 0) || !std::isfinite(h)) throw std::invalid_argument("SCF averaging width must be finite and positive");
	std::array<std::vector<std::pair<int, Real>>, 3> weights;
	for (int d = 0; d < 3; ++d) {
		Real const x = Real(physical[d] / length_) + (d == 0 ? diagnostics_.rotationCenter : 0);
		Real const lo = std::max(x-h/2, lower_), hi = std::min(x+h/2, lower_+n_*width_);
		if (hi <= lo) return {};
		int const begin = std::max(0, int(std::floor((lo-lower_)/width_)));
		int const end = std::min(n_-1, int(std::floor((hi-lower_)/width_)));
		for (int i = begin; i <= end; ++i) {
			Real const overlap = std::max(Real(0), std::min(hi, lower_+(i+1)*width_) - std::max(lo,lower_+i*width_));
			weights[d].emplace_back(i,overlap/h);
		}
	}
	Cell result;
	for (auto [k,wz] : weights[2]) for (auto [j,wy] : weights[1]) for (auto [i,wx] : weights[0]) {
		auto const& a = cells_[index(i,j,k,n_)]; Real const w = wx*wy*wz;
		result.density += w*a.density; result.pressure += w*a.pressure;
		for (int s = 0; s < 4; ++s) result.partial[s] += w*a.partial[s];
	}
	return result;
}

void BinaryScf::writeJson(std::ostream& out) const {
	auto const& d = diagnostics_;
	auto const& p = parameters_;
	Real const massUnit = units::value(p.primaryMass)*(1+p.massRatio);
	Real const length = units::value(length_);
	Real const energyUnit = units::value(constants::G)*massUnit*massUnit/length;
	Real const angularUnit = massUnit*length*length*units::value(units::sqrt(constants::G*density_));
	out << std::setprecision(17) << "{\n  \"cells_per_axis\": " << n_ << ",\n  \"iterations\": " << d.iterations
		<< ",\n  \"reference_width_in_separations\": " << p.referenceWidth
		<< ",\n  \"density_tolerance\": " << p.tolerance << ",\n  \"relaxation\": " << p.relaxation << ",\n  \"history\": " << p.history
		<< ",\n  \"hydro_gamma\": " << Real(5)/3 << ",\n  \"length_unit_cm\": " << length
		<< ",\n  \"density_unit_g_cm3\": " << units::value(density_) << ",\n  \"atmosphere_density_g_cm3\": " << atmosphere_*units::value(density_)
		<< ",\n  \"atmosphere_pressure_dyn_cm2\": " << atmosphere_*units::value(pressure_)
		<< ",\n  \"density_residual\": " << d.densityResidual << ",\n  \"bernoulli_residual\": " << d.bernoulliResidual
		<< ",\n  \"max_star_bernoulli_residual\": " << d.maxStarBernoulliResidual
		<< ",\n  \"virial_residual\": " << d.virialResidual << ",\n  \"virial_tolerance\": " << p.virialTolerance
		<< ",\n  \"separation_cm\": " << units::value(p.separation) << ",\n  \"omega_rad_s\": " << units::value(angularVelocity())
		<< ",\n  \"l1_x_cm\": " << (d.l1-d.rotationCenter)*length
		<< ",\n  \"kinetic_energy_erg\": " << d.kinetic*energyUnit << ",\n  \"potential_energy_erg\": " << d.potential*energyUnit
		<< ",\n  \"pressure_integral_erg\": " << d.pressureIntegral*energyUnit
		<< ",\n  \"orbital_angular_momentum_g_cm2_s\": " << d.orbitalAngularMomentum*angularUnit
		<< ",\n  \"spin_angular_momentum_g_cm2_s\": " << d.spinAngularMomentum*angularUnit << ",\n  \"stars\": [\n";
	for (int s = 0; s < 2; ++s) {
		out << "    {\"mass_g\": " << d.mass[s]*massUnit << ", \"core_mass_fraction\": " << d.coreMass[s]/d.mass[s]
			<< ", \"volume_equivalent_diameter_cm\": " << 2*std::cbrt(3*d.volume[s]/(4*std::numbers::pi))*length
			<< ", \"center_x_cm\": " << (d.center[s]-d.rotationCenter)*length << ", \"central_density_g_cm3\": " << d.centralDensity[s]*units::value(density_)
			<< ", \"core_index\": " << p.coreIndex[s] << ", \"envelope_index\": " << p.envelopeIndex[s]
			<< ", \"interface_density_fraction\": " << p.interfaceFraction[s] << ", \"density_jump\": " << p.densityJump[s]
			<< ", \"fill\": " << p.fill[s] << ", \"bernoulli_cm2_s2\": " << d.bernoulli[s]*energyUnit/massUnit << '}' << (s == 0 ? ",\n" : "\n");
	}
	out << "  ]\n}";
}
} // namespace octotigerII::problems
