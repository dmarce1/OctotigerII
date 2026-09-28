#pragma once
#include <array>
#include <cmath>
#include <octotigerII/helmholtz/helmholtz.hpp>
#include <stdexcept>
#include <tuple>
#include <vector>
namespace octotigerII::helmholtz::detail {
inline void validateState(double rho, double t, double a, double z) {
    if (!(std::isfinite(rho) && rho > 0 && std::isfinite(t) && t > 0 && std::isfinite(a) && a > 0 &&
          std::isfinite(z) && z > 0 && z <= a))
        throw std::invalid_argument("Helmholtz requires finite rho,T,Abar > 0 and 0 < Zbar <= Abar");
}
void validateGrid(const Grid &);
// One-based indexing intentionally follows the original Fortran.
struct Field {
    int ni{};
    std::vector<double> data;
    void resize(int n, int m) {
        ni = n + 1;
        data.resize((n + 1) * std::size_t(m + 1));
    }
    double &operator()(int i, int j) { return data[std::size_t(j) * ni + i]; }
    double operator()(int i, int j) const { return data[std::size_t(j) * ni + i]; }
};
struct Table {
    Grid grid;
    Field f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt;
    Field dpdf, dpdfd, dpdft, dpdfdt, ef, efd, eft, efdt, xf, xfd, xft, xfdt;
    std::vector<double> d, t, dt_sav, dt2_sav, dti_sav, dt2i_sav, dt3i_sav, dd_sav, dd2_sav, ddi_sav,
        dd2i_sav, dd3i_sav;
    explicit Table(Grid);
    std::array<Field *, 21> fields();
    auto bounds() const {
        return std::make_tuple(
            grid.densityPoints, grid.temperaturePoints, grid.logTemperatureMin, grid.logTemperatureMax,
            (grid.temperaturePoints - 1) / (grid.logTemperatureMax - grid.logTemperatureMin),
            grid.logDensityMin, grid.logDensityMax,
            (grid.densityPoints - 1) / (grid.logDensityMax - grid.logDensityMin));
    }
};
} // namespace octotigerII::helmholtz::detail
