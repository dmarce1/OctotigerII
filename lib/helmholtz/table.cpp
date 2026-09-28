// Table I/O follows F. X. Timmes's read_helm_table; see original below.
// Generation uses his direct EOS and the Timmes & Swesty free-energy method.
#include "direct.hpp"
#include "internal.hpp"
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
namespace octotigerII::helmholtz::detail {
void validateGrid(const Grid &g) {
    if (g.densityPoints < 2 || g.temperaturePoints < 2 || !std::isfinite(g.logDensityMin) ||
        !std::isfinite(g.logDensityMax) || !std::isfinite(g.logTemperatureMin) ||
        !std::isfinite(g.logTemperatureMax) || g.logDensityMin >= g.logDensityMax ||
        g.logTemperatureMin >= g.logTemperatureMax || g.logDensityMin < -12 || g.logDensityMax > 15 ||
        g.logTemperatureMin < 3 || g.logTemperatureMax > 13 ||
        std::size_t(g.densityPoints) * g.temperaturePoints > 10000000)
        throw std::invalid_argument("Invalid Helmholtz grid (supported logD [-12,15], logT [3,13])");
}
std::array<Field *, 21> Table::fields() {
    return {&f,     &fd,     &ft, &fdd, &ftt, &fdt,  &fddt, &fdtt, &fddtt, &dpdf, &dpdfd,
            &dpdft, &dpdfdt, &ef, &efd, &eft, &efdt, &xf,   &xfd,  &xft,   &xfdt};
}
Table::Table(Grid g) : grid(g) {
    validateGrid(g);
    for (auto *a : fields())
        a->resize(g.densityPoints, g.temperaturePoints);
    const auto axis = [](int n, double lo, double hi, std::vector<double> &x, std::vector<double> &h,
                         std::vector<double> &h2, std::vector<double> &inv, std::vector<double> &inv2,
                         std::vector<double> &inv3) {
        for (auto *a : {&x, &h, &h2, &inv, &inv2, &inv3})
            a->resize(n + 1);
        const double step = (hi - lo) / (n - 1);
        for (int i = 1; i <= n; ++i)
            x[i] = std::pow(10.0, lo + (i - 1) * step);
        for (int i = 1; i < n; ++i) {
            h[i] = x[i + 1] - x[i];
            if (!(h[i] > 0))
                throw std::invalid_argument("Helmholtz grid nodes are not distinct in double precision");
            h2[i] = h[i] * h[i];
            inv[i] = 1 / h[i];
            inv2[i] = 1 / h2[i];
            inv3[i] = inv2[i] * inv[i];
        }
    };
    axis(g.densityPoints, g.logDensityMin, g.logDensityMax, d, dd_sav, dd2_sav, ddi_sav, dd2i_sav, dd3i_sav);
    axis(g.temperaturePoints, g.logTemperatureMin, g.logTemperatureMax, t, dt_sav, dt2_sav, dti_sav, dt2i_sav,
         dt3i_sav);
}
} // namespace octotigerII::helmholtz::detail
namespace octotigerII::helmholtz {
Eos::Eos() : Eos(OCTOTIGERII_HELMHOLTZ_TABLE, Grid{}) {}
Eos::Eos(const std::string &path, Grid grid) {
    auto table = std::make_shared<detail::Table>(grid);
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open Helmholtz table: " + path);
    const auto fields = table->fields();
    for (auto block : {std::pair{0, 9}, std::pair{9, 13}, std::pair{13, 17}, std::pair{17, 21}})
        for (int j = 1; j <= grid.temperaturePoints; ++j)
            for (int i = 1; i <= grid.densityPoints; ++i)
                for (int k = block.first; k < block.second; ++k) {
                    // Accept both E and Fortran D exponents. Require the entire token.
                    std::string token;
                    if (!(in >> token))
                        throw std::runtime_error("Truncated Helmholtz table: " + path);
                    for (char &c : token)
                        if (c == 'D' || c == 'd')
                            c = 'e';
                    std::size_t used{};
                    const double v = std::stod(token, &used);
                    if (used != token.size() || !std::isfinite(v))
                        throw std::runtime_error("Invalid Helmholtz table value: " + path);
                    (*fields[k])(i, j) = v;
                }
    std::string extra;
    if (in >> extra)
        throw std::runtime_error("Extra data in Helmholtz table; check grid dimensions: " + path);
    table_ = std::move(table);
}
GenerationReport generateTable(const std::string &path, Grid grid, double step) {
    if (!std::isfinite(step) || step < 1e-5 || step > 1e-2)
        throw std::invalid_argument("Mixed-derivative step must be in [1e-5,1e-2]");
    detail::Table table(grid);
    GenerationReport report;
    const auto fields = table.fields();
    for (int j = 1; j <= grid.temperaturePoints; ++j) {
        const double t = table.t[j];
        for (int i = 1; i <= grid.densityPoints; ++i) {
            const double d = table.d[i], d2 = d * d, d3 = d2 * d;
            const auto q = detail::direct(d, t);
            if (q.st <= 0)
                ++report.nonpositiveElectronHeatCapacityNodes;
            if (q.pd < 0)
                ++report.negativeElectronCompressibilityNodes;
            // F=e-Ts; F_D=P/D^2; F_T=-s. Derivatives are at fixed D or T.
            table.f(i, j) = q.e - t * q.s;
            table.fd(i, j) = q.p / d2;
            table.ft(i, j) = -q.s;
            table.fdd(i, j) = q.pd / d2 - 2 * q.p / d3;
            table.ftt(i, j) = -q.st;
            table.fdt(i, j) = q.pt / d2;
            table.fddt(i, j) = q.pdt / d2 - 2 * q.pt / d3;
            table.fdtt(i, j) = q.ptt / d2;
            // The distributed direct EOS exposes second thermodynamic derivatives,
            // not P_DDT or F_DDTT. Obtain these remaining columns by a fourth-order
            // centered temperature derivative. This stencil is OctoII glue, not
            // a claim to reproduce Timmes's unpublished table-building driver.
            const auto m2 = detail::direct(d, t * (1 - 2 * step));
            const auto m1 = detail::direct(d, t * (1 - step));
            const auto p1 = detail::direct(d, t * (1 + step));
            const auto p2 = detail::direct(d, t * (1 + 2 * step));
            const auto fddt = [&](const detail::DirectPoint &a) { return a.pdt / d2 - 2 * a.pt / d3; };
            table.fddtt(i, j) = (fddt(m2) - 8 * fddt(m1) + 8 * fddt(p1) - fddt(p2)) / (12 * step * t);
            table.dpdf(i, j) = q.pd;
            table.dpdfd(i, j) = q.pdd;
            table.dpdft(i, j) = q.pdt;
            table.dpdfdt(i, j) = (m2.pdd - 8 * m1.pdd + 8 * p1.pdd - p2.pdd) / (12 * step * t);
            table.ef(i, j) = q.eta;
            table.efd(i, j) = q.etad;
            table.eft(i, j) = q.etat;
            table.efdt(i, j) = q.etadt;
            table.xf(i, j) = q.n;
            table.xfd(i, j) = q.nd;
            table.xft(i, j) = q.nt;
            table.xfdt(i, j) = q.ndt;
            for (auto *f : fields)
                if (!std::isfinite((*f)(i, j)))
                    throw std::runtime_error("Nonfinite generated table entry at D=" + std::to_string(d) +
                                             ", T=" + std::to_string(t));
        }
    }
    // Complete calculation before opening the output, so a solver failure
    // cannot replace an existing table with a partial numerical result.
    std::ofstream out(path);
    out << std::scientific << std::setprecision(17);
    for (auto block : {std::pair{0, 9}, std::pair{9, 13}, std::pair{13, 17}, std::pair{17, 21}})
        for (int j = 1; j <= grid.temperaturePoints; ++j)
            for (int i = 1; i <= grid.densityPoints; ++i) {
                for (int k = block.first; k < block.second; ++k)
                    out << (*fields[k])(i, j) << ' ';
                out << '\n';
            }
    out.close();
    if (!out)
        throw std::runtime_error("Failed writing Helmholtz table: " + path);
    std::ofstream metadata(path + ".meta.json");
    metadata << std::setprecision(17) << "{\n  \"format\": \"Timmes four-block ASCII\",\n"
             << "  \"density_points\": " << grid.densityPoints << ",\n"
             << "  \"temperature_points\": " << grid.temperaturePoints << ",\n"
             << "  \"log_density_min\": " << grid.logDensityMin << ",\n"
             << "  \"log_density_max\": " << grid.logDensityMax << ",\n"
             << "  \"log_temperature_min\": " << grid.logTemperatureMin << ",\n"
             << "  \"log_temperature_max\": " << grid.logTemperatureMax << ",\n"
             << "  \"mixed_derivative_step\": " << step << ",\n"
             << "  \"quadrature\": \"Timmes/Aparicio 20-point Legendre and Laguerre, IEEE double\",\n"
             << "  \"constants\": \"Timmes const.dek, 2006 CODATA\",\n"
             << "  \"source_sha256\": \"4c3e45924c00cff751885377c936cdd44793838c3b38db6ba4c5fec1c58710d3\",\n"
             << "  \"nonpositive_electron_heat_capacity_nodes\": "
             << report.nonpositiveElectronHeatCapacityNodes << ",\n"
             << "  \"negative_electron_compressibility_nodes\": "
             << report.negativeElectronCompressibilityNodes << ",\n"
             << "  \"accuracy_certified\": false\n}\n";
    metadata.close();
    if (!metadata)
        throw std::runtime_error("Failed writing Helmholtz metadata: " + path);
    return report;
}
} // namespace octotigerII::helmholtz

// Original read_helm_table, F. X. Timmes:
//      subroutine read_helm_table
//      include 'implno.dek'
//      include 'helm_table_storage.dek'
//
//! this routine reads the helmholtz eos file, and
//! must be called once before the helmeos routine is invoked.
//
//! declare local variables
//      integer          i,j
//      double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
//                       dd,dd2,ddi,dd2i,dd3i
//
//
//! open the file (use softlinks to input the desired table)
//
//       open(unit=19,file='helm_table.dat',status='old')
//
//
//! for standard table limits
//       tlo   = 3.0d0
//       thi   = 13.0d0
//       tstp  = (thi - tlo)/float(jmax-1)
//       tstpi = 1.0d0/tstp
//       dlo   = -12.0d0
//       dhi   = 15.0d0
//       dstp  = (dhi - dlo)/float(imax-1)
//       dstpi = 1.0d0/dstp
//
//! read the helmholtz free energy and its derivatives
//       do j=1,jmax
//        tsav = tlo + (j-1)*tstp
//        t(j) = 10.0d0**(tsav)
//        do i=1,imax
//         dsav = dlo + (i-1)*dstp
//         d(i) = 10.0d0**(dsav)
//         read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
//                  fddt(i,j),fdtt(i,j),fddtt(i,j)
//        enddo
//       enddo
//!       write(6,*) 'read main table'
//
//
//! read the pressure derivative with density table
//       do j=1,jmax
//        do i=1,imax
//         read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
//        enddo
//       enddo
//!       write(6,*) 'read dpdd table'
//
//! read the electron chemical potential table
//       do j=1,jmax
//        do i=1,imax
//         read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
//        enddo
//       enddo
//!       write(6,*) 'read eta table'
//
//! read the number density table
//       do j=1,jmax
//        do i=1,imax
//         read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
//        enddo
//       enddo
//!       write(6,*) 'read xne table'
//
//! close the file
//      close(unit=19)
//
//
//! construct the temperature and density deltas and their inverses
//       do j=1,jmax-1
//        dth          = t(j+1) - t(j)
//        dt2         = dth * dth
//        dti         = 1.0d0/dth
//        dt2i        = 1.0d0/dt2
//        dt3i        = dt2i*dti
//        dt_sav(j)   = dth
//        dt2_sav(j)  = dt2
//        dti_sav(j)  = dti
//        dt2i_sav(j) = dt2i
//        dt3i_sav(j) = dt3i
//       end do
//       do i=1,imax-1
//        dd          = d(i+1) - d(i)
//        dd2         = dd * dd
//        ddi         = 1.0d0/dd
//        dd2i        = 1.0d0/dd2
//        dd3i        = dd2i*ddi
//        dd_sav(i)   = dd
//        dd2_sav(i)  = dd2
//        ddi_sav(i)  = ddi
//        dd2i_sav(i) = dd2i
//        dd3i_sav(i) = dd3i
//       enddo
//
//
//
//!      write(6,*)
//!      write(6,*) 'finished reading eos table'
//!      write(6,04) 'imax=',imax,' jmax=',jmax
//!04    format(1x,4(a,i4))
//!      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
//!      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
//!03    format(1x,4(a,1pe11.3))
//!      write(6,*)
//
//      return
//      end
//
//
//
//
//
//
