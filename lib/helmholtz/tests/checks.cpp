#include "constants.hpp"
#include "direct.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <octotigerII/helmholtz/helmholtz.hpp>
#include <stdexcept>
using namespace octotigerII::helmholtz;
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
void close(double x, double y, double tol, const char *message) {
    require(std::isfinite(x) && std::isfinite(y) && std::abs(x - y) <= tol * std::max(std::abs(y), 1e-200),
            message);
}
template <class F> void rejects(F f) {
    bool caught = false;
    try {
        f();
    } catch (const std::exception &) {
        caught = true;
    }
    require(caught, "Expected rejection");
}
int main() {
    const auto path = std::filesystem::current_path() / "helmholtz-test-table.dat";
    try {
        // Default data selection must work independently of the working directory.
        const Eos supplied;
        const auto cold = supplied.evaluate(2e9, 1e4, 12, 6);
        require(std::isfinite(cold.ptot) && cold.ptot > 0 && std::isfinite(cold.deept) && cold.cv > 0,
                "Supplied table cold-degenerate state");
        // Classical electron gas: independently known P=nkT and e=3kT/(2mu).
        const auto classical = detail::direct(1e-6, 1e6);
        close(classical.p, 1e-6 * detail::avo * detail::kerg * 1e6, 5e-4, "Classical pressure");
        close(classical.e, 1.5 * detail::avo * detail::kerg * 1e6, 5e-4, "Classical energy");
        // Nonrelativistic zero-temperature electron pressure, with a small
        // finite-temperature/relativistic correction allowed at this state.
        const double rho = 1e3, n = rho * detail::avo;
        const double p0 = std::pow(3 * detail::pi * detail::pi, 2.0 / 3) * detail::hbar * detail::hbar /
                          (5 * detail::me) * std::pow(n, 5.0 / 3);
        close(detail::direct(rho, 1e4).p, p0, 0.01, "Degenerate pressure");
        Grid g{21, 21, 4, 5, 7, 8};
        const auto report = generateTable(path.string(), g);
        require(report.nonpositiveElectronHeatCapacityNodes == 0 &&
                    report.negativeElectronCompressibilityNodes == 0,
                "Generated table sign diagnostics");
        Eos eos(path.string(), g);
        for (double d : {1e4, 1.618e4, 4.72e4, 1e5}) {
            for (double t : {1e7, 1.337e7, 6.82e7, 1e8}) {
                const auto q = eos.evaluate(2 * d, t, 12, 6);
                const auto reference = detail::direct(d, t);
                close(q.pele, reference.p, 2e-6, "Generated table pressure");
                close(q.eele, 0.5 * reference.e, 2e-6, "Generated table energy");
                close(q.sele, 0.5 * reference.s, 2e-6, "Generated table entropy");
                close(q.dpepd, 0.5 * reference.pd, 2e-5, "Generated pressure derivative");
                require(std::abs(q.dse) < 2e-12 && std::abs(q.dpe) < 2e-12 && std::abs(q.dsp) < 2e-12,
                        "Maxwell relations");
                require(q.det > 0 && q.gam1 > 0 && q.cs > 0, "Stable thermodynamics");
                close(q.ptot, q.pgas + q.prad, 1e-14, "Photon decomposition");
            }
        }
        // Compare an actual centered perturbation to the returned derivative.
        const double rho1 = 60000, t1 = 2.718e7, h = 1e-5;
        auto q = eos.evaluate(rho1, t1, 12, 6);
        close((eos.evaluate(rho1, t1 * (1 + h), 12, 6).ptot - eos.evaluate(rho1, t1 * (1 - h), 12, 6).ptot) /
                  (2 * h * t1),
              q.dpt, 1e-6, "Temperature derivative");
        // Halving the mixed-derivative stencil step should leave the
        // thermodynamics unchanged well below the interpolation error.
        generateTable(path.string(), g, 5e-4);
        Eos halfStep(path.string(), g);
        const auto half = halfStep.evaluate(rho1, t1, 12, 6);
        close(half.pele, q.pele, 1e-10, "Mixed derivative step: pressure");
        close(half.eele, q.eele, 1e-10, "Mixed derivative step: energy");
        close(half.deept, q.deept, 1e-8, "Mixed derivative step: heat capacity");
        // Resolve the pair-production transition with successively finer tables.
        const double pairD = 6131.90350185317, pairT = 1653076520.037604;
        const auto pair = detail::direct(pairD, pairT);
        double previous = 1.0;
        for (int intervals : {10, 20, 40}) {
            Grid pairGrid{2 * intervals + 1, intervals + 1, 3, 5, 9, 10};
            generateTable(path.string(), pairGrid);
            Eos pairEos(path.string(), pairGrid);
            const auto value = pairEos.evaluate(2 * pairD, pairT, 12, 6);
            const double error = std::abs(2 * value.dpepd / pair.pd - 1);
            std::cout << "Pair compressibility error, " << intervals << " intervals/decade: " << error
                      << '\n';
            require(error < previous / 4, "Pair compressibility table convergence");
            previous = error;
        }
        rejects([&] { eos.evaluate(0, t1, 12, 6); });
        rejects([&] { eos.evaluate(rho1, std::numeric_limits<double>::quiet_NaN(), 12, 6); });
        rejects([&] { eos.evaluate(rho1, t1, 0, 6); });
        rejects([&] { eos.evaluate(rho1, t1, 12, 13); });
        rejects([&] { eos.evaluate(1, t1, 12, 6); });
        rejects([&] { eos.evaluate(rho1, 1e9, 12, 6); });
        rejects([&] { Eos wrong(path.string(), {2, 2, 4, 5, 7, 8}); });
        {
            std::ofstream truncated(path);
            truncated << "1 2 3\n";
        }
        rejects([&] { Eos bad(path.string(), g); });
        std::filesystem::remove(path);
        std::filesystem::remove(path.string() + ".meta.json");
        std::cout << "Helmholtz analytic, generated-table, derivative, consistency and input checks passed\n";
    } catch (const std::exception &e) {
        std::filesystem::remove(path);
        std::filesystem::remove(path.string() + ".meta.json");
        std::cerr << e.what() << '\n';
        return 1;
    }
}
