#include "octotigerII/hydro/helmholtzClosure.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <octotigerII/helmholtz/helmholtz.hpp>
#include <stdexcept>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/synchronization/mutex.hpp>
#endif

namespace octotigerII::hydro {
namespace {
    // Immutable tables shared per locality, including short-lived HydroSystem objects.
    // Initialization is serialized with an HPX-aware mutex in HPX builds.
    std::shared_ptr<const helmholtz::Eos> table(std::string const& path) {
#ifdef OCTOTIGERII_WITH_HPX
        static hpx::mutex mutex;
#else
        static std::mutex mutex;
#endif
        static std::map<std::string, std::shared_ptr<const helmholtz::Eos>> tables;
        std::lock_guard lock(mutex);
        auto& result = tables[path];
        if (!result) result = path.empty() ? std::make_shared<helmholtz::Eos>() : std::make_shared<helmholtz::Eos>(path);
        return result;
    }
}    // namespace
HelmholtzClosure::HelmholtzClosure(std::string const& path, bool gasOnly, double floor)
  : eos_(table(path))
  , gasOnly_(gasOnly)
  , floor_(floor) {
    if (!std::isfinite(floor) || floor < 1e3 || floor >= 1e13) throw std::invalid_argument("Helmholtz temperature floor must be in [1000, 1e13) K");
}
HelmholtzClosure::Point HelmholtzClosure::at(double rho, double temperature, double abar, double zbar) const {
    auto const q = eos_->evaluate(rho, temperature, abar, zbar);
    return {temperature, gasOnly_ ? q.egas : q.etot, gasOnly_ ? q.pgas : q.ptot, gasOnly_ ? q.sgas : q.stot, gasOnly_ ? q.degast : q.det,
        gasOnly_ ? q.dpgast : q.dpt, (gasOnly_ ? q.gam1_gas * q.pgas : q.gam1 * q.ptot) / rho};
}
HelmholtzClosure::Point HelmholtzClosure::invert(double rho, double target, double abar, double zbar, Variable variable, bool floor) const {
    if (!std::isfinite(target)) throw std::invalid_argument("Nonfinite Helmholtz inversion target");
    auto value = [&](Point const& q) { return variable == Variable::Energy ? q.energy : variable == Variable::Pressure ? q.pressure : q.entropy; };
    auto derivative = [&](Point const& q) { return variable == Variable::Energy ? q.cv : variable == Variable::Pressure ? q.dpdt : q.cv / q.temperature; };
    auto low = at(rho, floor_, abar, zbar), high = at(rho, 1e13, abar, zbar);
    if (target <= value(low)) {
        if (!floor && target < value(low) * (1 - 1e-12)) throw std::out_of_range("Helmholtz target below temperature floor");
        return low;
    }
    if (target > value(high)) throw std::out_of_range("Helmholtz target above maximum table temperature");
    double t = std::sqrt(low.temperature * high.temperature);
    for (int iteration = 0; iteration < 100; ++iteration) {
        auto const q = at(rho, t, abar, zbar);
        double const residual = value(q) - target;
        if (residual == 0 || high.temperature / low.temperature - 1 < 4e-15) return q;
        if (residual < 0)
            low = q;
        else
            high = q;
        double const slope = derivative(q);
        double const next = slope > 0 ? t - residual / slope : -1;
        // Keep a sizeable bracket contraction even in the degenerate limit.
        double const lo = std::exp(0.9 * std::log(low.temperature) + 0.1 * std::log(high.temperature));
        double const hi = std::exp(0.1 * std::log(low.temperature) + 0.9 * std::log(high.temperature));
        t = slope > 0 && std::isfinite(next) && next > lo && next < hi ? next : std::sqrt(low.temperature * high.temperature);
    }
    throw std::runtime_error("Helmholtz temperature inversion did not converge");
}
}    // namespace octotigerII::hydro
