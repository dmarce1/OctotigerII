#include "octotigerII/gravity/ewald.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#include "octotigerII/profiling.hpp"

namespace octotigerII::gravity::ewald {
namespace {
	using Scalar = long double;
	using Vector = std::array<Scalar, 3>;
	constexpr Scalar pi = 3.141592653589793238462643383279502884L;
	int index(int x, int y, int z) {
		int const n = x + y + z;
		return n * n + (z == 0 ? x : n + 1 + x);
	}
	Scalar power(Scalar x, int n) {
		Scalar result = 1;
		while (n--)
			result *= x;
		return result;
	}
	std::vector<diagonal::Index> indices(int p) {
		std::vector<diagonal::Index> result;
		for (int n = 0; n <= p; ++n) {
			for (int x = 0; x <= n; ++x)
				result.push_back({x, n - x, 0});
			for (int x = 0; x < n; ++x)
				result.push_back({x, n - 1 - x, 1});
		}
		return result;
	}
	std::vector<std::pair<Scalar, Scalar>> quadrature(int n, Scalar alpha) {
		using std::abs;
		using std::cos;
		std::vector<std::pair<Scalar, Scalar>> result;
		for (int i = 0; i < n; ++i) {
			Scalar z = cos(pi * (i + 0.75L) / (n + 0.5L)), dp = 0;
			for (int it = 0; it < 80; ++it) {
				Scalar p = 1, previous = 0;
				for (int j = 1; j <= n; ++j) {
					Scalar const old = previous;
					previous = p;
					p = ((2 * j - 1) * z * previous - (j - 1) * old) / j;
				}
				dp = n * (z * p - previous) / (z * z - 1);
				Scalar const dz = p / dp;
				z -= dz;
				if (abs(dz) < 4 * std::numeric_limits<Scalar>::epsilon()) break;
				if (it == 79) throw std::runtime_error("Ewald Gauss-Legendre quadrature failed");
			}
			result.emplace_back(alpha * (1 + z) / 2, alpha / ((1 - z * z) * dp * dp));
		}
		return result;
	}

	// d^n exp(-a x^2)/dx^n, generated without numerical differentiation.
	std::vector<Scalar> gaussian(Scalar x, Scalar a, int degree) {
		using std::exp;
		std::vector<Scalar> g(degree + 1);
		g[0] = exp(-a * x * x);
		if (degree) g[1] = -2 * a * x * g[0];
		for (int n = 1; n < degree; ++n)
			g[n + 1] = -2 * a * (x * g[n] + n * g[n - 1]);
		return g;
	}
	Scalar phase(Scalar cosine, Scalar sine, int degree) {
		switch (degree % 4) {
		case 0:
			return cosine;
		case 1:
			return -sine;
		case 2:
			return -cosine;
		default:
			return sine;
		}
	}

	// Cartesian derivatives of a radial function with D_{k+1}=-D'_k/r.
	void radial(std::vector<Scalar>& out, Vector r, std::vector<Scalar> const& D, std::vector<diagonal::Index> const& ix, Scalar sign) {
		for (std::size_t slot = 0; slot < ix.size(); ++slot) {
			auto const a = ix[slot];
			Scalar sum = 0, cx = 1;
			for (int i = 0; i <= a.x / 2; ++i) {
				Scalar cy = 1;
				for (int j = 0; j <= a.y / 2; ++j) {
					int const k = a.degree() - i - j;
					sum += ((k & 1) ? -1 : 1) * cx * cy * power(r[0], a.x - 2 * i) * power(r[1], a.y - 2 * j) * power(r[2], a.z) * D[k];
					cy *= Scalar((a.y - 2 * j) * (a.y - 2 * j - 1)) / (2 * (j + 1));
				}
				cx *= Scalar((a.x - 2 * i) * (a.x - 2 * i - 1)) / (2 * (i + 1));
			}
			out[slot] += sign * sum;
		}
	}

	std::map<std::array<int, 7>, std::shared_ptr<Operator const>> cache;
	std::map<std::array<int, 6>, std::array<double, 4>> pointCache;
	std::shared_mutex cacheMutex, pointMutex;
}	 // namespace

double Derivatives::at(int x, int y, int z) const {
	if (z >= 2) return -at(x + 2, y, z - 2) - at(x, y + 2, z - 2) + ((x == 0 && y == 0 && z == 2) ? laplacian : 0);
	return values.at(index(x, y, z));
}

Derivatives derivatives(diagonal::Vector input, diagonal::Vector lengths, int degree, double alphaInput, int points) {
	using std::abs;
	using std::ceil;
	using std::cos;
	using std::erf;
	using std::erfc;
	using std::exp;
	using std::expm1;
	using std::isfinite;
	using std::log;
	using std::sin;
	using std::sqrt;
	if (degree < 0 || degree > 20 || points < 16 || points > 256 || !isfinite(alphaInput) || alphaInput < 0)
		throw std::invalid_argument("Invalid Ewald derivative controls");
	int dimensions = 0, freeAxis = -1, periodicAxis = -1;
	Scalar minimum = std::numeric_limits<Scalar>::infinity(), volume = 1;
	Vector r{}, period{};
	for (int d = 0; d < 3; ++d) {
		if (!isfinite(input[d]) || !isfinite(lengths[d]) || lengths[d] < 0) throw std::invalid_argument("Invalid Ewald geometry");
		r[d] = input[d];
		period[d] = lengths[d];
		if (period[d]) {
			++dimensions;
			volume *= period[d];
			minimum = std::min(minimum, period[d]);
			periodicAxis = d;
		} else
			freeAxis = d;
	}
	if (!dimensions) throw std::invalid_argument("Ewald needs at least one periodic axis");
	Scalar const alpha = alphaInput > 0 ? Scalar(alphaInput) : 2 / minimum;
	Scalar const alpha2 = alpha * alpha;
	auto const ix = indices(degree);
	std::vector<Scalar> out(ix.size(), 0), D(degree + 1, 0);
	auto const quad = quadrature(points, alpha);
	Scalar const r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
	// Add erf(alpha*r)/r directly, regular at the physical self point.
	for (auto [t, w] : quad) {
		Scalar value = 2 / sqrt(pi) * w * exp(-t * t * r2);
		for (int k = 0; k <= degree; ++k) {
			D[k] += value;
			value *= 2 * t * t;
		}
	}
	radial(out, r, D, ix, 1);
	// Exponentially small truncation; cutoff increases with derivative order.
	Scalar const cutoff = 8 + Scalar(degree) / 8;
	std::array<int, 3> bound{};
	for (int d = 0; d < 3; ++d)
		if (period[d]) bound[d] = int(ceil((cutoff / alpha + abs(r[d])) / period[d]));
	for (int x = -bound[0]; x <= bound[0]; ++x)
		for (int y = -bound[1]; y <= bound[1]; ++y)
			for (int z = -bound[2]; z <= bound[2]; ++z) {
				if (!x && !y && !z) continue;
				Vector const v{r[0] - x * period[0], r[1] - y * period[1], r[2] - z * period[2]};
				Scalar const s2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
				if (s2 == 0) throw std::invalid_argument("Ewald correction evaluated at an excluded image singularity");
				if (alpha2 * s2 > cutoff * cutoff) continue;
				Scalar const s = sqrt(s2);
				D[0] = erfc(alpha * s) / s;
				Scalar gaussianTerm = 2 * alpha / sqrt(pi) * exp(-alpha2 * s2);
				for (int k = 0; k < degree; ++k) {
					D[k + 1] = ((2 * k + 1) * D[k] + gaussianTerm) / s2;
					gaussianTerm *= 2 * alpha2;
				}
				radial(out, v, D, ix, -1);
			}

	// Nonzero reciprocal modes. Opposite modes are combined into one cosine.
	for (int d = 0; d < 3; ++d)
		bound[d] = period[d] ? int(ceil(alpha * period[d] * cutoff / pi)) : 0;
	for (int x = -bound[0]; x <= bound[0]; ++x)
		for (int y = -bound[1]; y <= bound[1]; ++y)
			for (int z = -bound[2]; z <= bound[2]; ++z) {
				if (x < 0 || (!x && y < 0) || (!x && !y && z <= 0)) continue;
				Vector k{period[0] ? 2 * pi * x / period[0] : 0, period[1] ? 2 * pi * y / period[1] : 0, period[2] ? 2 * pi * z / period[2] : 0};
				Scalar const k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
				if (k2 / (4 * alpha2) > cutoff * cutoff) continue;
				Scalar const angle = k[0] * r[0] + k[1] * r[1] + k[2] * r[2], cs = cos(angle), sn = sin(angle);
				std::array<std::vector<Scalar>, 3> powers;
				for (int d = 0; d < 3; ++d) {
					powers[d].resize(degree + 1, 1);
					for (int n = 1; n <= degree; ++n)
						powers[d][n] = powers[d][n - 1] * k[d];
				}
				if (dimensions == 3) {
					Scalar const weight = -8 * pi / volume / k2 * exp(-k2 / (4 * alpha2));
					for (std::size_t i = 0; i < ix.size(); ++i) {
						auto a = ix[i];
						out[i] += weight * powers[0][a.x] * powers[1][a.y] * powers[2][a.z] * phase(cs, sn, a.degree());
					}
				} else if (dimensions == 2) {
					Scalar const km = sqrt(k2), u = r[freeAxis];
					Scalar plus = exp(km * u) * erfc(km / (2 * alpha) + alpha * u);
					Scalar minus = exp(-km * u) * erfc(km / (2 * alpha) - alpha * u);
					auto q = gaussian(u, alpha2, degree);
					std::vector<Scalar> H(degree + 1);
					Scalar const factor = 2 * alpha / sqrt(pi) * exp(-k2 / (4 * alpha2));
					for (int n = 0; n <= degree; ++n) {
						H[n] = plus + minus;
						plus = km * plus - factor * q[n];
						minus = -km * minus + factor * q[n];
					}
					for (std::size_t i = 0; i < ix.size(); ++i) {
						auto a = ix[i];
						std::array<int, 3> e{a.x, a.y, a.z};
						int const f = e[freeAxis];
						e[freeAxis] = 0;
						out[i] += -2 * pi / volume / km * H[f] * powers[0][e[0]] * powers[1][e[1]] * powers[2][e[2]] * phase(cs, sn, a.degree() - f);
					}
				} else {
					for (auto [t, w] : quad) {
						std::array<std::vector<Scalar>, 3> g;
						for (int d = 0; d < 3; ++d)
							g[d] = period[d] ? powers[d] : gaussian(r[d], t * t, degree);
						Scalar const weight = -4 / volume * w / t * exp(-k2 / (4 * t * t));
						for (std::size_t i = 0; i < ix.size(); ++i) {
							auto a = ix[i];
							std::array<int, 3> e{a.x, a.y, a.z};
							out[i] += weight * g[0][a.x] * g[1][a.y] * g[2][a.z] * phase(cs, sn, e[periodicAxis]);
						}
					}
				}
			}
	// Retain the transverse zero mode for 1P and 2P (line/sheet gravity).
	if (dimensions == 3)
		out[0] += pi / (volume * alpha2);
	else if (dimensions == 2) {
		Scalar const u = r[freeAxis];
		auto g = gaussian(u, alpha2, degree);
		for (std::size_t i = 0; i < ix.size(); ++i) {
			auto a = ix[i];
			std::array<int, 3> e{a.x, a.y, a.z};
			int const f = e[freeAxis];
			if (f != a.degree()) continue;
			if (f == 0)
				out[i] += 2 * sqrt(pi) / volume * (g[0] / alpha + sqrt(pi) * u * erf(alpha * u));
			else if (f == 1)
				out[i] += 2 * pi / volume * erf(alpha * u);
			else
				out[i] += 4 * sqrt(pi) * alpha / volume * g[f - 2];
		}
	} else {
		// Finite-part gauge log(alpha*minimum). Its alpha derivative cancels
		// that of the real sum; no source/background density is subtracted.
		out[0] -= 2 / volume * log(alpha * minimum);
		for (auto [t, w] : quad) {
			std::array<std::vector<Scalar>, 3> g;
			Scalar transverse2 = 0;
			for (int d = 0; d < 3; ++d) {
				g[d] = period[d] ? std::vector<Scalar>(degree + 1, 1) : gaussian(r[d], t * t, degree);
				if (!period[d]) transverse2 += r[d] * r[d];
			}
			for (std::size_t i = 0; i < ix.size(); ++i) {
				auto a = ix[i];
				std::array<int, 3> e{a.x, a.y, a.z};
				if (e[periodicAxis]) continue;
				Scalar const value = i == 0 ? expm1(-t * t * transverse2) : g[0][a.x] * g[1][a.y] * g[2][a.z];
				out[i] -= 2 / volume * w / t * value;
			}
		}
	}
	Derivatives result;
	result.laplacian = dimensions == 3 ? double(-4 * pi / volume) : 0;
	for (auto value : out) {
		if (!isfinite(double(value))) throw std::runtime_error("Nonfinite Ewald derivative");
		result.values.push_back(double(value));
	}
	return result;
}

Operator::Operator(int p, diagonal::Offset r, diagonal::Offset periods)
  : count_(diagonal::coefficientCount(p)) {
	using std::pow;
	profiling::Region profile("gravity.ewald.operator_setup");
	double scale = std::numeric_limits<double>::infinity();
	for (int v : periods)
		if (v > 0) scale = std::min(scale, double(v));
	diagonal::Vector x{}, lengths{};
	for (int d = 0; d < 3; ++d) {
		x[d] = r[d] / scale;
		lengths[d] = periods[d] / scale;
	}
	auto const kernel = derivatives(x, lengths, 2 * p);
	laplace_ = kernel.laplacian / pow(scale, 3);
	matrix_.resize(std::size_t(count_) * count_);
	auto const& ix = diagonal::indices(p);
	for (int i = 0; i < count_; ++i)
		for (int j = 0; j < count_; ++j) {
			auto const a = ix[i], b = ix[j];
			matrix_[std::size_t(i) * count_ + j] =
				((b.degree() & 1) ? -1 : 1) * kernel.at(a.x + b.x, a.y + b.y, a.z + b.z) / pow(scale, 1 + a.degree() + b.degree());
		}
}
void Operator::add(diagonal::Coefficients& l, diagonal::Coefficients const& m, double h) const {
	if (l.size() != std::size_t(count_ + 1) || (m.size() != 1 && m.size() != std::size_t(count_ + 1)) || !(h > 0))
		throw std::invalid_argument("Invalid Ewald M2L data");
	std::size_t const count = m.size() == 1 ? 1 : std::size_t(count_);
	for (int i = 0; i < count_; ++i) {
		double sum = 0;
		for (std::size_t j = 0; j < count; ++j)
			sum += matrix_[std::size_t(i) * count_ + j] * m[j];
		l[i] += sum / h;
	}
	if (m.size() > 1) l[0] += laplace_ * m.back() / h;
	l.back() += laplace_ * m[0] / h;
}
std::shared_ptr<Operator const> getOperator(int p, diagonal::Offset r, diagonal::Offset periods) {
	std::array<int, 7> const key{p, r[0], r[1], r[2], periods[0], periods[1], periods[2]};
	{
		std::shared_lock lock(cacheMutex);
		auto it = cache.find(key);
		if (it != cache.end()) return it->second;
	}
	// Expensive setup must not hold a lock needed by other worker tasks.
	auto op = std::make_shared<Operator const>(p, r, periods);
	std::unique_lock lock(cacheMutex);
	return cache.emplace(key, std::move(op)).first->second;
}
void addDirect(diagonal::Coefficients& l, double mass, diagonal::Offset r, diagonal::Offset periods, double h) {
	using std::pow;
	if (mass == 0) return;
	std::array<int, 6> const key{r[0], r[1], r[2], periods[0], periods[1], periods[2]};
	std::array<double, 4> value{};
	bool found = false;
	{
		std::shared_lock lock(pointMutex);
		auto it = pointCache.find(key);
		if (it != pointCache.end()) {
			value = it->second;
			found = true;
		}
	}
	if (!found) {
		double scale = std::numeric_limits<double>::infinity();
		for (int v : periods)
			if (v > 0) scale = std::min(scale, double(v));
		diagonal::Vector x{}, lengths{};
		for (int d = 0; d < 3; ++d) {
			x[d] = r[d] / scale;
			lengths[d] = periods[d] / scale;
		}
		auto kernel = derivatives(x, lengths, 1);
		for (int i = 0; i < 4; ++i)
			value[i] = kernel.values[i] / pow(scale, i == 0 ? 1 : 2);
		std::unique_lock lock(pointMutex);
		pointCache.emplace(key, value);
	}
	for (int i = 0; i < 4; ++i)
		l[i] += mass / h * value[i];
}
}	 // namespace octotigerII::gravity::ewald
