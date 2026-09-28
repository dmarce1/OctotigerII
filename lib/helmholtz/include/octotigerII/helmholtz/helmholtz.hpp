// C++ interface to F. X. Timmes's Helmholtz EOS. See NOTICE.md.
#pragma once
#include <memory>
#include <string>
namespace octotigerII::helmholtz {
namespace detail {
struct Table;
}
// Logarithmic grid in D = rho * Zbar/Abar [g/cm^3] and T [K].
// Defaults match the user-supplied Timmes table (20 intervals per decade).
struct Grid {
    int densityPoints = 541;
    int temperaturePoints = 201;
    double logDensityMin = -12;
    double logDensityMax = 15;
    double logTemperatureMin = 3;
    double logTemperatureMax = 13;
};
// Timmes's output names, with the Fortran _row suffix removed.
// CGS: p [erg/cm^3], e [erg/g], s [erg/g/K], cs [cm/s].
// Derivative suffixes t,d,a,z denote T,rho,Abar,Zbar at fixed other inputs.
// `ele` contains the combined electron/positron contribution; `pos` is zero
// in the tabulated evaluator, exactly as in the supplied Fortran.
// `gas` excludes photons. cs is Timmes's relativistic sound speed;
// Newtonian hydro should use sqrt(gam1*p/rho) for its chosen components.
struct Result {
    double ptot{};
    double dpt{};
    double dpd{};
    double dpa{};
    double dpz{};
    double etot{};
    double det{};
    double ded{};
    double dea{};
    double dez{};
    double stot{};
    double dst{};
    double dsd{};
    double dsa{};
    double dsz{};
    double pgas{};
    double dpgast{};
    double dpgasd{};
    double dpgasa{};
    double dpgasz{};
    double egas{};
    double degast{};
    double degasd{};
    double degasa{};
    double degasz{};
    double sgas{};
    double dsgast{};
    double dsgasd{};
    double dsgasa{};
    double dsgasz{};
    double prad{};
    double dpradt{};
    double dpradd{};
    double dprada{};
    double dpradz{};
    double erad{};
    double deradt{};
    double deradd{};
    double derada{};
    double deradz{};
    double srad{};
    double dsradt{};
    double dsradd{};
    double dsrada{};
    double dsradz{};
    double pion{};
    double dpiont{};
    double dpiond{};
    double dpiona{};
    double dpionz{};
    double eion{};
    double deiont{};
    double deiond{};
    double deiona{};
    double deionz{};
    double sion{};
    double dsiont{};
    double dsiond{};
    double dsiona{};
    double dsionz{};
    double xni{};
    double pele{};
    double ppos{};
    double dpept{};
    double dpepd{};
    double dpepa{};
    double dpepz{};
    double eele{};
    double epos{};
    double deept{};
    double deepd{};
    double deepa{};
    double deepz{};
    double sele{};
    double spos{};
    double dsept{};
    double dsepd{};
    double dsepa{};
    double dsepz{};
    double xnem{};
    double xne{};
    double dxnet{};
    double dxned{};
    double dxnea{};
    double dxnez{};
    double xnp{};
    double zeff{};
    double etaele{};
    double detat{};
    double detad{};
    double detaa{};
    double detaz{};
    double etapos{};
    double pcou{};
    double dpcout{};
    double dpcoud{};
    double dpcoua{};
    double dpcouz{};
    double ecou{};
    double decout{};
    double decoud{};
    double decoua{};
    double decouz{};
    double scou{};
    double dscout{};
    double dscoud{};
    double dscoua{};
    double dscouz{};
    double plasg{};
    double dse{};
    double dpe{};
    double dsp{};
    double cv_gas{};
    double cp_gas{};
    double gam1_gas{};
    double gam2_gas{};
    double gam3_gas{};
    double nabad_gas{};
    double cs_gas{};
    double cv{};
    double cp{};
    double gam1{};
    double gam2{};
    double gam3{};
    double nabad{};
    double cs{};
};
class Eos {
  public:
    // Explicit initialization; no lazy I/O, shared mutable state, or locks.
    // Tables use Timmes's four-block ASCII format, with an explicit grid.
    // Load the checksum-verified Timmes table bundled with this source build.
    Eos();
    explicit Eos(const std::string &tablePath, Grid grid = {});
    Result evaluate(double density, double temperature, double abar, double zbar) const;

  private:
    std::shared_ptr<const detail::Table> table_;
};
// Build the table from direct Fermi-Dirac quadrature, never from a seed table.
// This is offline work; the evaluator needs only the resulting file.
// mixedDerivativeStep is the relative temperature step for the fourth
// free-energy derivative and third pressure derivative (5-point stencil).
struct GenerationReport {
    // Sign checks are diagnostics, not a complete accuracy certification.
    int nonpositiveElectronHeatCapacityNodes{};
    int negativeElectronCompressibilityNodes{};
};
// Also writes PATH.meta.json with the grid, numerical method and sign checks.
// The original direct quadrature can lose thermal accuracy in cold degenerate
// states. Inspect this report and validate the intended physical domain.
GenerationReport generateTable(const std::string &path, Grid grid = {}, double mixedDerivativeStep = 1e-3);
} // namespace octotigerII::helmholtz
