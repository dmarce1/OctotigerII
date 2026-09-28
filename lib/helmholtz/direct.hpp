#pragma once
namespace octotigerII::helmholtz::detail {
// Fully ionized electron/positron gas, rho here means rho*Ye (Ye=1).
struct DirectPoint {
    double p, e, s, pd, pt, pdd, pdt, ptt, st, sdd;
    double eta, etad, etat, etadt, n, nd, nt, ndt;
};
DirectPoint direct(double density, double temperature);
} // namespace octotigerII::helmholtz::detail
