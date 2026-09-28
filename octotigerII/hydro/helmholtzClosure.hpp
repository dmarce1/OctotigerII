// C++ hydrodynamic closure for F. X. Timmes's EOS; upstream credit in lib/helmholtz/NOTICE.md.
#pragma once
namespace octotigerII::helmholtz {
class Eos;
}
#include <memory>
#include <string>

namespace octotigerII::hydro {
class HelmholtzClosure {
public:
    struct Point
    {
        double temperature, energy, pressure, entropy, cv, dpdt, soundSquared;
    };
    enum class Variable
    {
        Energy,
        Pressure,
        Entropy
    };
    HelmholtzClosure(std::string const& path, bool gasOnly, double floor);
    Point at(double rho, double temperature, double abar, double zbar) const;
    Point invert(double rho, double target, double abar, double zbar, Variable variable, bool floor = true) const;
    bool gasOnly() const {
        return gasOnly_;
    }
    double temperatureFloor() const {
        return floor_;
    }

private:
    std::shared_ptr<const helmholtz::Eos> eos_;
    bool gasOnly_;
    double floor_;
};
}    // namespace octotigerII::hydro
