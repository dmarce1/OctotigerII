// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/gravity/diagonal/fmm.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#include <tuple>

namespace octotigerII::gravity::diagonal {
namespace {
constexpr double pi = 3.1415926535897932384626433832795;
int index(int x, int y, int z) {
	const int n = x + y + z;
	return n * n + (z == 0 ? x : n + 1 + x);
}
double power(double x, int n) {
	double r = 1;
	for (int i = 0; i < n; ++i)
		r *= x;
	return r;
}
double factorial(int n) {
	double r = 1;
	for (int i = 2; i <= n; ++i)
		r *= i;
	return r;
}
void reduce(Coefficients& a, int x, int y, int z, double v) {
	if (z < 2)
		a[index(x, y, z)] += v;
	else {
		reduce(a, x + 2, y, z - 2, -v);
		reduce(a, x, y + 2, z - 2, -v);
	}
}
std::map<std::array<int, 4>, std::shared_ptr<const Operator>> cache;
std::shared_mutex cacheMutex;
// Golub-Welsch is unnecessary at these tiny orders: Newton iteration on L_n.
std::vector<std::pair<double, double>> laguerre(int n) {
	std::vector<std::pair<double, double>> q;
	double z = 0;
	for (int i = 0; i < n; ++i) {
		if (i == 0)
			z = 3.0 / (1 + 2.4 * n);
		else if (i == 1)
			z += 15.0 / (1 + 2.5 * n);
		else {
			const double a = i - 1;
			z += (1 + 2.55 * a) / (1.9 * a) * (z - q[i - 2].first);
		}
		double prev = 0, deriv = 0;
		bool converged = false;
		for (int it = 0; it < 50; ++it) {
			double l = 1;
			prev = 0;
			for (int j = 1; j <= n; ++j) {
				const double old = prev;
				prev = l;
				l = ((2 * j - 1 - z) * prev - (j - 1) * old) / j;
			}
			deriv = n * (l - prev) / z;
			const double dz = l / deriv;
			z -= dz;
			if (std::abs(dz) < 2e-15 * (1 + std::abs(z))) {
				converged = true;
				break;
			}
		}
		if (!converged)
			throw std::runtime_error("diagonal FMM Laguerre quadrature failed");
		// Reevaluate derivative at the converged root.
		double l = 1;
		prev = 0;
		for (int j = 1; j <= n; ++j) {
			const double old = prev;
			prev = l;
			l = ((2 * j - 1 - z) * prev - (j - 1) * old) / j;
		}
		deriv = n * (l - prev) / z;
		q.emplace_back(z, 1 / (z * deriv * deriv));
	}
	return q;
}
} // namespace
int coefficientCount(int p) {
	if (p < 0 || p > 5)
		throw std::invalid_argument("diagonal FMM order must be 0..5");
	return (p + 1) * (p + 1);
}
const std::vector<Index>& indices(int p) {
	coefficientCount(p);
	static const auto tables = [] {
		std::array<std::vector<Index>, 6> result;
		for (int p = 0; p <= 5; ++p)
			for (int n = 0; n <= p; ++n) {
				for (int x = 0; x <= n; ++x)
					result[p].push_back({x, n - x, 0});
				for (int x = 0; x < n; ++x)
					result[p].push_back({x, n - 1 - x, 1});
			}
		return result;
	}();
	return tables[p];
}
Coefficients point(double mass, const Vector& s, int p) {
	Coefficients a(coefficientCount(p), 0);
	for (int n = 0; n <= p; ++n)
		for (int x = 0; x <= n; ++x)
			for (int y = 0; y <= n - x; ++y) {
				const int z = n - x - y;
				reduce(a, x, y, z,
					   mass * power(s[0], x) * power(s[1], y) * power(s[2], z) /
						   (factorial(x) * factorial(y) * factorial(z)));
			}
	return a;
}
Coefficients shiftMultipole(const Coefficients& a, const Vector& s, double ratio, int p) {
	const auto e = point(1, s, p);
	const auto& ix = indices(p);
	if (a.size() != 1 && a.size() != ix.size())
		throw std::invalid_argument("invalid multipole size");
	Coefficients b(ix.size(), 0);
	for (std::size_t i = 0; i < a.size(); ++i)
		for (std::size_t j = 0; j < ix.size(); ++j) {
			if (ix[i].degree() + ix[j].degree() > p)
				continue;
			reduce(b, ix[i].x + ix[j].x, ix[i].y + ix[j].y, ix[i].z + ix[j].z,
				   a[i] * power(ratio, ix[i].degree()) * e[j]);
		}
	return b;
}
double derivative(const Coefficients& l, int x, int y, int z) {
	if (z >= 2)
		return -derivative(l, x + 2, y, z - 2) - derivative(l, x, y + 2, z - 2);
	return l.at(index(x, y, z));
}
Coefficients shiftLocal(const Coefficients& l, const Vector& s, double ratio, int p) {
	const auto& ix = indices(p);
	if (l.size() != ix.size())
		throw std::invalid_argument("invalid local size");
	Coefficients b(ix.size(), 0);
	for (std::size_t i = 0; i < ix.size(); ++i) {
		const auto a = ix[i];
		for (int n = 0; n <= p - a.degree(); ++n)
			for (int x = 0; x <= n; ++x)
				for (int y = 0; y <= n - x; ++y) {
					const int z = n - x - y;
					b[i] += derivative(l, a.x + x, a.y + y, a.z + z) * power(s[0], x) *
							power(s[1], y) * power(s[2], z) /
							(factorial(x) * factorial(y) * factorial(z));
				}
		b[i] *= power(ratio, a.degree());
	}
	return b;
}
std::array<double, 4> evaluate(const Coefficients& l, const Vector& s, int p) {
	const auto b = shiftLocal(l, s, 1, p);
	return {b[0], derivative(b, 1, 0, 0), derivative(b, 0, 1, 0), derivative(b, 0, 0, 1)};
}
Operator::Operator(int p, Offset r) : count_(coefficientCount(p)) {
	if (p < 3)
		throw std::invalid_argument("diagonal FMM M2L order must be 3,4,5");
	const double radius = std::hypot(double(r[0]), double(r[1]), double(r[2]));
	if (radius == 0)
		throw std::invalid_argument("zero M2L separation");
	Vector e{r[0] / radius, r[1] / radius, r[2] / radius};
	int a = 0;
	for (int d = 1; d < 3; ++d)
		if (std::abs(e[d]) < std::abs(e[a]))
			a = d;
	Vector u{};
	u[a] = 1;
	const double dot = e[a];
	for (int d = 0; d < 3; ++d)
		u[d] -= dot * e[d];
	const double norm = std::hypot(u[0], u[1], u[2]);
	for (auto& v : u)
		v /= norm;
	Vector v{e[1] * u[2] - e[2] * u[1], e[2] * u[0] - e[0] * u[2], e[0] * u[1] - e[1] * u[0]};
	// Source degree p + local derivative degree p: radial degree <=2p,
	// angular Fourier degree <=2p. Both quadratures are exact at these sizes.
	const int angles = 2 * p + 1;
	for (auto [t, w] : laguerre(p + 1))
		for (int j = 0; j < angles; ++j) {
			const double angle = 2 * pi * j / angles;
			std::array<std::complex<double>, 3> k;
			for (int d = 0; d < 3; ++d)
				k[d] = t / radius *
					   std::complex<double>(-e[d], u[d] * std::cos(angle) + v[d] * std::sin(angle));
			// exp(k.R)=exp(-t) is the diagonal translation. Its exp(-t)
			// factor is already in the Gauss-Laguerre weight w.
			diagonal_.push_back(-w / (radius * angles));
			for (const auto b : indices(p)) {
				std::complex<double> z = 1;
				for (int n = 0; n < b.x; ++n)
					z *= k[0];
				for (int n = 0; n < b.y; ++n)
					z *= k[1];
				for (int n = 0; n < b.z; ++n)
					z *= k[2];
				fromWave_.push_back(z);
				toWave_.push_back((b.degree() % 2 ? -1.0 : 1.0) * z);
			}
		}
}
void Operator::add(Coefficients& l, const Coefficients& m, double h) const {
	if (!(h > 0) || !std::isfinite(h) || (m.size() != 1 && m.size() != std::size_t(count_)) ||
		(l.size() != 4 && l.size() != std::size_t(count_)))
		throw std::invalid_argument("invalid M2L data");
	for (std::size_t q = 0; q < diagonal_.size(); ++q) {
		std::complex<double> wave = 0;
		const auto off = q * count_;
		for (std::size_t i = 0; i < m.size(); ++i)
			wave += toWave_[off + i] * m[i];
		wave *= diagonal_[q] / h;
		for (std::size_t i = 0; i < l.size(); ++i)
			l[i] += (fromWave_[off + i] * wave).real();
	}
}
std::shared_ptr<const Operator> getOperator(int p, Offset r) {
	const std::array<int, 4> key{p, r[0], r[1], r[2]};
	{
		std::shared_lock lock(cacheMutex);
		const auto it = cache.find(key);
		if (it != cache.end())
			return it->second;
	}
	// Short construction, no suspension or HPX calls under the lock. Build once.
	std::unique_lock lock(cacheMutex);
	auto it = cache.find(key);
	if (it != cache.end())
		return it->second;
	auto op = std::make_shared<const Operator>(p, r);
	cache.emplace(key, op);
	return op;
}
std::size_t cachedOperatorCount() {
	std::shared_lock lock(cacheMutex);
	return cache.size();
}
void addDirect(Coefficients& l, double m, Offset r, double h) {
	const double distance = std::hypot(double(r[0]), double(r[1]), double(r[2]));
	if (distance == 0 || !(h > 0))
		throw std::invalid_argument("invalid P2P separation");
	l.at(0) -= m / (h * distance);
	const double f = m / (h * distance * distance * distance);
	l.at(index(1, 0, 0)) += f * r[0];
	l.at(index(0, 1, 0)) += f * r[1];
	l.at(index(0, 0, 1)) += f * r[2];
}
} // namespace octotigerII::gravity::diagonal
