#include "octotigerII/gravity/isolatedPotential.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace octotigerII::gravity {
namespace {
	using Complex = std::complex<Real>;
	std::size_t index(int x, int y, int z, int n) { return (std::size_t(z) * n + y) * n + x; }
	// Iterative radix-two transform. All lines use one precomputed twiddle table;
	// the inverse normalization is applied once per dimension.
	void transform(std::vector<Complex>& a, int n, bool inverse) {
		std::vector<Complex> roots(n/2), line(n);
		for (int i = 0; i < n/2; ++i) roots[i] = std::polar(Real(1), (inverse ? 2 : -2) * std::numbers::pi_v<Real> * i / n);
		for (int d = 0; d < 3; ++d) {
			std::size_t const stride = d == 0 ? 1 : d == 1 ? n : n*n;
			for (int u = 0; u < n; ++u) for (int v = 0; v < n; ++v) {
				std::size_t const base = d == 0 ? index(0,u,v,n) : d == 1 ? index(u,0,v,n) : index(u,v,0,n);
				for (int i = 0; i < n; ++i) line[i] = a[base + i*stride];
				for (int i = 1, j = 0; i < n; ++i) {
					int bit = n >> 1;
					for (; j & bit; bit >>= 1) j ^= bit;
					j ^= bit;
					if (i < j) std::swap(line[i], line[j]);
				}
				for (int length = 2; length <= n; length *= 2) {
					int const half = length/2;
					for (int first = 0; first < n; first += length) for (int j = 0; j < half; ++j) {
						auto const even = line[first+j], odd = roots[j*(n/length)] * line[first+j+half];
						line[first+j] = even+odd; line[first+j+half] = even-odd;
					}
				}
				for (int i = 0; i < n; ++i) a[base+i*stride] = inverse ? line[i]/Real(n) : line[i];
			}
		}
	}
}
IsolatedPotential::IsolatedPotential(int n, Real dx) : n_(n), padded_(0) {
	if (n < 2 || n > 128 || (n & (n-1)) || !std::isfinite(dx*dx) || !(dx > 0) || dx*dx == 0)
		throw std::invalid_argument("Isolated potential requires power-of-two cells=2..128 and finite positive cell width");
	padded_ = 2*n;
	kernel_.resize(std::size_t(padded_)*padded_*padded_);
	for (int z = 0; z < padded_; ++z) for (int y = 0; y < padded_; ++y) for (int x = 0; x < padded_; ++x) {
		Real const r = std::hypot(Real(std::min(x,padded_-x)), Real(std::min(y,padded_-y)), Real(std::min(z,padded_-z)));
		kernel_[index(x,y,z,padded_)] = r > 0 ? -dx*dx/r : 0;
	}
	transform(kernel_, padded_, false);
}
std::vector<Real> IsolatedPotential::operator()(std::vector<Real> const& density) const {
	if (density.size() != std::size_t(n_)*n_*n_) throw std::invalid_argument("Isolated potential density size mismatch");
	std::vector<Complex> work(kernel_.size());
	for (int z = 0; z < n_; ++z) for (int y = 0; y < n_; ++y) for (int x = 0; x < n_; ++x) {
		Real const rho = density[index(x,y,z,n_)];
		if (!std::isfinite(rho) || rho < 0) throw std::invalid_argument("Isolated potential requires finite nonnegative density");
		work[index(x,y,z,padded_)] = rho;
	}
	transform(work,padded_,false);
	for (std::size_t i = 0; i < work.size(); ++i) work[i] *= kernel_[i];
	transform(work,padded_,true);
	std::vector<Real> result(density.size());
	for (int z = 0; z < n_; ++z) for (int y = 0; y < n_; ++y) for (int x = 0; x < n_; ++x)
		result[index(x,y,z,n_)] = work[index(x,y,z,padded_)].real();
	return result;
}
} // namespace octotigerII::gravity
