// C++ translation of F. X. (Frank) Timmes's stellar EOS routines.
// Original scientific algorithms: Timmes & Arnett (1999), ApJS 125, 277;
// Timmes & Swesty (2000), ApJS 126, 501. See NOTICE.md for provenance.
// Original Fortran statements are retained immediately above their translations.
#include "constants.hpp"
#include "internal.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace octotigerII::helmholtz {
Result Eos::evaluate(double den, double temp, double abar, double zbar) const {
    using namespace detail;
    validateState(den, temp, abar, zbar);
    const auto &table = *table_;
    const auto [imax, jmax, tlo, thi, tstpi, dlo, dhi, dstpi] = table.bounds();
    Result out{};
    [[maybe_unused]] double ytot1{};
    [[maybe_unused]] double ye{};
    [[maybe_unused]] double x{};
    [[maybe_unused]] double y{};
    [[maybe_unused]] double zz{};
    [[maybe_unused]] double zzi{};
    [[maybe_unused]] double deni{};
    [[maybe_unused]] double tempi{};
    [[maybe_unused]] double xni{};
    [[maybe_unused]] double dxnidd{};
    [[maybe_unused]] double dxnida{};
    [[maybe_unused]] double dpepdt{};
    [[maybe_unused]] double dpepdd{};
    [[maybe_unused]] double deepdt{};
    [[maybe_unused]] double deepdd{};
    [[maybe_unused]] double dsepdd{};
    [[maybe_unused]] double dsepdt{};
    [[maybe_unused]] double dpraddd{};
    [[maybe_unused]] double dpraddt{};
    [[maybe_unused]] double deraddd{};
    [[maybe_unused]] double deraddt{};
    [[maybe_unused]] double dpiondd{};
    [[maybe_unused]] double dpiondt{};
    [[maybe_unused]] double deiondd{};
    [[maybe_unused]] double deiondt{};
    [[maybe_unused]] double dsraddd{};
    [[maybe_unused]] double dsraddt{};
    [[maybe_unused]] double dsiondd{};
    [[maybe_unused]] double dsiondt{};
    [[maybe_unused]] double dse{};
    [[maybe_unused]] double dpe{};
    [[maybe_unused]] double dsp{};
    [[maybe_unused]] double kt{};
    [[maybe_unused]] double ktinv{};
    [[maybe_unused]] double prad{};
    [[maybe_unused]] double erad{};
    [[maybe_unused]] double srad{};
    [[maybe_unused]] double pion{};
    [[maybe_unused]] double eion{};
    [[maybe_unused]] double sion{};
    [[maybe_unused]] double xnem{};
    [[maybe_unused]] double pele{};
    [[maybe_unused]] double eele{};
    [[maybe_unused]] double sele{};
    [[maybe_unused]] double pres{};
    [[maybe_unused]] double ener{};
    [[maybe_unused]] double entr{};
    [[maybe_unused]] double dpresdd{};
    [[maybe_unused]] double dpresdt{};
    [[maybe_unused]] double denerdd{};
    [[maybe_unused]] double denerdt{};
    [[maybe_unused]] double dentrdd{};
    [[maybe_unused]] double dentrdt{};
    [[maybe_unused]] double cv{};
    [[maybe_unused]] double cp{};
    [[maybe_unused]] double gam1{};
    [[maybe_unused]] double gam2{};
    [[maybe_unused]] double gam3{};
    [[maybe_unused]] double chit{};
    [[maybe_unused]] double chid{};
    [[maybe_unused]] double nabad{};
    [[maybe_unused]] double sound{};
    [[maybe_unused]] double etaele{};
    [[maybe_unused]] double detadt{};
    [[maybe_unused]] double detadd{};
    [[maybe_unused]] double xnefer{};
    [[maybe_unused]] double dxnedt{};
    [[maybe_unused]] double dxnedd{};
    [[maybe_unused]] double s{};
    [[maybe_unused]] double pgas{};
    [[maybe_unused]] double dpgasdd{};
    [[maybe_unused]] double dpgasdt{};
    [[maybe_unused]] double dpgasda{};
    [[maybe_unused]] double dpgasdz{};
    [[maybe_unused]] double egas{};
    [[maybe_unused]] double degasdd{};
    [[maybe_unused]] double degasdt{};
    [[maybe_unused]] double degasda{};
    [[maybe_unused]] double degasdz{};
    [[maybe_unused]] double sgas{};
    [[maybe_unused]] double dsgasdd{};
    [[maybe_unused]] double dsgasdt{};
    [[maybe_unused]] double dsgasda{};
    [[maybe_unused]] double dsgasdz{};
    [[maybe_unused]] double cv_gas{};
    [[maybe_unused]] double cp_gas{};
    [[maybe_unused]] double gam1_gas{};
    [[maybe_unused]] double gam2_gas{};
    [[maybe_unused]] double gam3_gas{};
    [[maybe_unused]] double chit_gas{};
    [[maybe_unused]] double chid_gas{};
    [[maybe_unused]] double nabad_gas{};
    [[maybe_unused]] double sound_gas{};
    [[maybe_unused]] double dpradda{};
    [[maybe_unused]] double deradda{};
    [[maybe_unused]] double dsradda{};
    [[maybe_unused]] double dpionda{};
    [[maybe_unused]] double deionda{};
    [[maybe_unused]] double dsionda{};
    [[maybe_unused]] double dpepda{};
    [[maybe_unused]] double deepda{};
    [[maybe_unused]] double dsepda{};
    [[maybe_unused]] double dpresda{};
    [[maybe_unused]] double denerda{};
    [[maybe_unused]] double dentrda{};
    [[maybe_unused]] double detada{};
    [[maybe_unused]] double dxneda{};
    [[maybe_unused]] double dpraddz{};
    [[maybe_unused]] double deraddz{};
    [[maybe_unused]] double dsraddz{};
    [[maybe_unused]] double dpiondz{};
    [[maybe_unused]] double deiondz{};
    [[maybe_unused]] double dsiondz{};
    [[maybe_unused]] double dpepdz{};
    [[maybe_unused]] double deepdz{};
    [[maybe_unused]] double dsepdz{};
    [[maybe_unused]] double dpresdz{};
    [[maybe_unused]] double denerdz{};
    [[maybe_unused]] double dentrdz{};
    [[maybe_unused]] double detadz{};
    [[maybe_unused]] double dxnedz{};
    [[maybe_unused]] int iat{};
    [[maybe_unused]] int jat{};
    [[maybe_unused]] double free{};
    [[maybe_unused]] double df_d{};
    [[maybe_unused]] double df_t{};
    [[maybe_unused]] double df_dd{};
    [[maybe_unused]] double df_tt{};
    [[maybe_unused]] double df_dt{};
    [[maybe_unused]] double xt{};
    [[maybe_unused]] double xd{};
    [[maybe_unused]] double mxt{};
    [[maybe_unused]] double mxd{};
    [[maybe_unused]] double si0t{};
    [[maybe_unused]] double si1t{};
    [[maybe_unused]] double si2t{};
    [[maybe_unused]] double si0mt{};
    [[maybe_unused]] double si1mt{};
    [[maybe_unused]] double si2mt{};
    [[maybe_unused]] double si0d{};
    [[maybe_unused]] double si1d{};
    [[maybe_unused]] double si2d{};
    [[maybe_unused]] double si0md{};
    [[maybe_unused]] double si1md{};
    [[maybe_unused]] double si2md{};
    [[maybe_unused]] double dsi0t{};
    [[maybe_unused]] double dsi1t{};
    [[maybe_unused]] double dsi2t{};
    [[maybe_unused]] double dsi0mt{};
    [[maybe_unused]] double dsi1mt{};
    [[maybe_unused]] double dsi2mt{};
    [[maybe_unused]] double dsi0d{};
    [[maybe_unused]] double dsi1d{};
    [[maybe_unused]] double dsi2d{};
    [[maybe_unused]] double dsi0md{};
    [[maybe_unused]] double dsi1md{};
    [[maybe_unused]] double dsi2md{};
    [[maybe_unused]] double ddsi0t{};
    [[maybe_unused]] double ddsi1t{};
    [[maybe_unused]] double ddsi2t{};
    [[maybe_unused]] double ddsi0mt{};
    [[maybe_unused]] double ddsi1mt{};
    [[maybe_unused]] double ddsi2mt{};
    [[maybe_unused]] double ddsi0d{};
    [[maybe_unused]] double ddsi1d{};
    [[maybe_unused]] double ddsi2d{};
    [[maybe_unused]] double ddsi0md{};
    [[maybe_unused]] double ddsi1md{};
    [[maybe_unused]] double ddsi2md{};
    [[maybe_unused]] double din{};
    double fi[37]{};
    [[maybe_unused]] double dsdd{};
    [[maybe_unused]] double dsda{};
    [[maybe_unused]] double lami{};
    [[maybe_unused]] double inv_lami{};
    [[maybe_unused]] double lamida{};
    [[maybe_unused]] double lamidd{};
    [[maybe_unused]] double plasg{};
    [[maybe_unused]] double plasgdd{};
    [[maybe_unused]] double plasgdt{};
    [[maybe_unused]] double plasgda{};
    [[maybe_unused]] double plasgdz{};
    [[maybe_unused]] double ecoul{};
    [[maybe_unused]] double decouldd{};
    [[maybe_unused]] double decouldt{};
    [[maybe_unused]] double decoulda{};
    [[maybe_unused]] double decouldz{};
    [[maybe_unused]] double pcoul{};
    [[maybe_unused]] double dpcouldd{};
    [[maybe_unused]] double dpcouldt{};
    [[maybe_unused]] double dpcoulda{};
    [[maybe_unused]] double dpcouldz{};
    [[maybe_unused]] double scoul{};
    [[maybe_unused]] double dscouldd{};
    [[maybe_unused]] double dscouldt{};
    [[maybe_unused]] double dscoulda{};
    [[maybe_unused]] double dscouldz{};
    double z{};
    const double sioncon = ((((2.0 * pi) * amu) * kerg) / (h * h));
    const double forth = (4.0 / 3.0);
    [[maybe_unused]] const double forpi = (4.0 * pi);
    const double kergavo = (kerg * avo);
    [[maybe_unused]] const double ikavo = (1.0 / kergavo);
    const double asoli3 = (asol / 3.0);
    const double light2 = (clight * clight);
    const double a1 = (-0.898004);
    const double b1 = 0.96786;
    const double c1 = 0.220703;
    const double d1 = (-0.86097);
    const double e1 = 2.5269;
    const double a2 = 0.29561;
    const double b2 = 1.9885;
    const double c2 = 0.288675;
    const double third = (1.0 / 3.0);
    const double esqu = (qe * qe);
    //      subroutine helmeos
    //      include 'implno.dek'
    //      include 'const.dek'
    //      include 'vector_eos.dek'
    //      include 'helm_table_storage.dek'
    //
    //
    //! given a temperature temp [K], density den [g/cm**3], and a composition
    //! characterized by abar and zbar, this routine returns most of the other
    //! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
    //! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
    //! their derivatives with respect to temperature, density, abar, and zbar.
    //! other quantites such the normalized chemical potential eta (plus its
    //! derivatives), number density of electrons and positron pair (along
    //! with their derivatives), adiabatic indices, specific heats, and
    //! relativistically correct sound speed are also returned.
    //!
    //! this routine assumes planckian photons, an ideal gas of ions,
    //! and an electron-positron gas with an arbitrary degree of relativity
    //! and degeneracy. interpolation in a table of the helmholtz free energy
    //! is used to return the electron-positron thermodynamic quantities.
    //! all other derivatives are analytic.
    //!
    //! references: cox & giuli chapter 24 ; timmes & swesty apj 1999
    //
    //
    //! declare
    //      integer          i,j
    //      double precision temp,den,abar,zbar,ytot1,ye, &
    //                       x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
    //                       dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
    //                       dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
    //                       deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
    //                       dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
    //                       sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
    //                       dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
    //                       gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
    //                       detadt,detadd,xnefer,dxnedt,dxnedd,s
    //
    //      double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
    //                       egas,degasdd,degasdt,degasda,degasdz, &
    //                       sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
    //                       cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
    //                       chit_gas,chid_gas,nabad_gas,sound_gas
    //
    //
    //      double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
    //      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
    //                        forth   = 4.0d0/3.0d0, &
    //                        forpi   = 4.0d0 * pi, &
    //                        kergavo = kerg * avo, &
    //                        ikavo   = 1.0d0/kergavo, &
    //                        asoli3  = asol/3.0d0, &
    //                        light2  = clight * clight)
    //
    //! for the abar derivatives
    //      double precision dpradda,deradda,dsradda, &
    //                       dpionda,deionda,dsionda, &
    //                       dpepda,deepda,dsepda, &
    //                       dpresda,denerda,dentrda, &
    //                       detada,dxneda
    //
    //! for the zbar derivatives
    //      double precision dpraddz,deraddz,dsraddz, &
    //                       dpiondz,deiondz,dsiondz, &
    //                       dpepdz,deepdz,dsepdz, &
    //                       dpresdz,denerdz,dentrdz, &
    //                       detadz,dxnedz
    //
    //! for the interpolations
    //      integer          iat,jat
    //      double precision free,df_d,df_t,df_dd,df_tt,df_dt
    //      double precision xt,xd,mxt,mxd, &
    //                       si0t,si1t,si2t,si0mt,si1mt,si2mt, &
    //                       si0d,si1d,si2d,si0md,si1md,si2md, &
    //                       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
    //                       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
    //                       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
    //                       ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
    //                       z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
    //                       dpsi2,ddpsi2,din,h5,fi(36), &
    //                       xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
    //                       w0t,w1t,w2t,w0mt,w1mt,w2mt, &
    //                       w0d,w1d,w2d,w0md,w1md,w2md
    //
    //
    //! for the uniform background coulomb correction
    //      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
    //                       plasg,plasgdd,plasgdt,plasgda,plasgdz, &
    //                       ecoul,decouldd,decouldt,decoulda,decouldz, &
    //                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
    //                       scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
    //                       a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
    //      parameter        (a1    = -0.898004d0, &
    //                        b1    =  0.96786d0, &
    //                        c1    =  0.220703d0, &
    //                        d1    = -0.86097d0, &
    //                        e1    =  2.5269d0, &
    //                        a2    =  0.29561d0, &
    //                        b2    =  1.9885d0, &
    //                        c2    =  0.288675d0, &
    //                        third =  1.0d0/3.0d0, &
    //                        esqu  =  qe * qe)
    //
    //
    //! quintic hermite polynomial statement functions
    //! psi0 and its derivatives
    //      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
    const auto psi0 = [&](double z) { return (((z * z * z) * ((z * (((-6.0) * z) + 15.0)) - 10.0)) + 1.0); };
    //      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
    const auto dpsi0 = [&](double z) { return ((z * z) * ((z * (((-30.0) * z) + 60.0)) - 30.0)); };
    //      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)
    const auto ddpsi0 = [&](double z) { return (z * ((z * (((-120.0) * z) + 180.0)) - 60.0)); };
    //
    //
    //! psi1 and its derivatives
    //      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
    const auto psi1 = [&](double z) { return (z * (((z * z) * ((z * (((-3.0) * z) + 8.0)) - 6.0)) + 1.0)); };
    //      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
    const auto dpsi1 = [&](double z) { return (((z * z) * ((z * (((-15.0) * z) + 32.0)) - 18.0)) + 1.0); };
    //      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
    const auto ddpsi1 = [&](double z) { return (z * ((z * (((-60.0) * z) + 96.0)) - 36.0)); };
    //
    //
    //! psi2  and its derivatives
    //      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
    const auto psi2 = [&](double z) { return (((0.5 * z) * z) * ((z * ((z * ((-z) + 3.0)) - 3.0)) + 1.0)); };
    //      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
    const auto dpsi2 = [&](double z) {
        return ((0.5 * z) * ((z * ((z * (((-5.0) * z) + 12.0)) - 9.0)) + 2.0));
    };
    //      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)
    const auto ddpsi2 = [&](double z) { return (0.5 * ((z * ((z * (((-20.0) * z) + 36.0)) - 18.0)) + 2.0)); };
    //
    //
    //! biquintic hermite polynomial statement function
    //      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
    //             fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
    //           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
    //           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
    //           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
    //           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
    //           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
    //           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
    //           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
    //           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
    //           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
    //           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
    //           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
    //           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
    //           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
    //           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
    //           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
    //           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
    //           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
    const auto h5 = [&](double w0t, double w1t, double w2t, double w0mt, double w1mt, double w2mt, double w0d,
                        double w1d, double w2d, double w0md, double w1md, double w2md) {
        return (((((((((((((((((((((((((((((((((((((fi[1] * w0d) * w0t) + ((fi[2] * w0md) * w0t)) +
                                                 ((fi[3] * w0d) * w0mt)) +
                                                ((fi[4] * w0md) * w0mt)) +
                                               ((fi[5] * w0d) * w1t)) +
                                              ((fi[6] * w0md) * w1t)) +
                                             ((fi[7] * w0d) * w1mt)) +
                                            ((fi[8] * w0md) * w1mt)) +
                                           ((fi[9] * w0d) * w2t)) +
                                          ((fi[10] * w0md) * w2t)) +
                                         ((fi[11] * w0d) * w2mt)) +
                                        ((fi[12] * w0md) * w2mt)) +
                                       ((fi[13] * w1d) * w0t)) +
                                      ((fi[14] * w1md) * w0t)) +
                                     ((fi[15] * w1d) * w0mt)) +
                                    ((fi[16] * w1md) * w0mt)) +
                                   ((fi[17] * w2d) * w0t)) +
                                  ((fi[18] * w2md) * w0t)) +
                                 ((fi[19] * w2d) * w0mt)) +
                                ((fi[20] * w2md) * w0mt)) +
                               ((fi[21] * w1d) * w1t)) +
                              ((fi[22] * w1md) * w1t)) +
                             ((fi[23] * w1d) * w1mt)) +
                            ((fi[24] * w1md) * w1mt)) +
                           ((fi[25] * w2d) * w1t)) +
                          ((fi[26] * w2md) * w1t)) +
                         ((fi[27] * w2d) * w1mt)) +
                        ((fi[28] * w2md) * w1mt)) +
                       ((fi[29] * w1d) * w2t)) +
                      ((fi[30] * w1md) * w2t)) +
                     ((fi[31] * w1d) * w2mt)) +
                    ((fi[32] * w1md) * w2mt)) +
                   ((fi[33] * w2d) * w2t)) +
                  ((fi[34] * w2md) * w2t)) +
                 ((fi[35] * w2d) * w2mt)) +
                ((fi[36] * w2md) * w2mt));
    };
    //
    //
    //
    //! cubic hermite polynomial statement functions
    //! psi0 & derivatives
    //      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
    const auto xpsi0 = [&](double z) { return (((z * z) * ((2.0 * z) - 3.0)) + 1.0); };
    //      xdpsi0(z) = z * (6.0d0*z - 6.0d0)
    const auto xdpsi0 = [&](double z) { return (z * ((6.0 * z) - 6.0)); };
    //
    //
    //! psi1 & derivatives
    //      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
    const auto xpsi1 = [&](double z) { return (z * ((z * (z - 2.0)) + 1.0)); };
    //      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0
    const auto xdpsi1 = [&](double z) { return ((z * ((3.0 * z) - 4.0)) + 1.0); };
    //
    //
    //! bicubic hermite polynomial statement function
    //      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
    //             fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
    //           + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
    //           + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
    //           + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
    //           + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
    //           + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
    //           + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
    //           + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt
    const auto h3 = [&](double w0t, double w1t, double w0mt, double w1mt, double w0d, double w1d, double w0md,
                        double w1md) {
        return (((((((((((((((((fi[1] * w0d) * w0t) + ((fi[2] * w0md) * w0t)) + ((fi[3] * w0d) * w0mt)) +
                            ((fi[4] * w0md) * w0mt)) +
                           ((fi[5] * w0d) * w1t)) +
                          ((fi[6] * w0md) * w1t)) +
                         ((fi[7] * w0d) * w1mt)) +
                        ((fi[8] * w0md) * w1mt)) +
                       ((fi[9] * w1d) * w0t)) +
                      ((fi[10] * w1md) * w0t)) +
                     ((fi[11] * w1d) * w0mt)) +
                    ((fi[12] * w1md) * w0mt)) +
                   ((fi[13] * w1d) * w1t)) +
                  ((fi[14] * w1md) * w1t)) +
                 ((fi[15] * w1d) * w1mt)) +
                ((fi[16] * w1md) * w1mt));
    };
    //
    //
    //
    //! popular format statements
    //01    format(1x,5(a,1pe11.3))
    //02    format(1x,a,1p4e16.8)
    //03    format(1x,4(a,1pe11.3))
    //04    format(1x,4(a,i4))
    //
    //
    //
    //! start of pipeline loop, normal execution starts here
    //      eosfail = .false.
    //      do j=jlo_eos,jhi_eos
    //
    //!       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
    //!       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'
    //
    //       temp  = temp_row(j)
    //       den   = den_row(j)
    //       abar  = abar_row(j)
    //       zbar  = zbar_row(j)
    //       ytot1 = 1.0d0/abar
    ytot1 = (1.0 / abar);
    //       ye    = max(1.0d-16,ytot1 * zbar)
    ye = std::max(1e-16, (ytot1 * zbar));
    //
    //
    //
    //! initialize
    //       deni    = 1.0d0/den
    deni = (1.0 / den);
    //       tempi   = 1.0d0/temp
    tempi = (1.0 / temp);
    //       kt      = kerg * temp
    kt = (kerg * temp);
    //       ktinv   = 1.0d0/kt
    ktinv = (1.0 / kt);
    //
    //
    //! radiation section:
    //       prad    = asoli3 * temp * temp * temp * temp
    prad = ((((asoli3 * temp) * temp) * temp) * temp);
    //       dpraddd = 0.0d0
    dpraddd = 0.0;
    //       dpraddt = 4.0d0 * prad*tempi
    dpraddt = ((4.0 * prad) * tempi);
    //       dpradda = 0.0d0
    dpradda = 0.0;
    //       dpraddz = 0.0d0
    dpraddz = 0.0;
    //
    //       erad    = 3.0d0 * prad*deni
    erad = ((3.0 * prad) * deni);
    //       deraddd = -erad*deni
    deraddd = ((-erad) * deni);
    //       deraddt = 3.0d0 * dpraddt*deni
    deraddt = ((3.0 * dpraddt) * deni);
    //       deradda = 0.0d0
    deradda = 0.0;
    //       deraddz = 0.0d0
    deraddz = 0.0;
    //
    //       srad    = (prad*deni + erad)*tempi
    srad = (((prad * deni) + erad) * tempi);
    //       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
    dsraddd = ((((dpraddd * deni) - ((prad * deni) * deni)) + deraddd) * tempi);
    //       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
    dsraddt = ((((dpraddt * deni) + deraddt) - srad) * tempi);
    //       dsradda = 0.0d0
    dsradda = 0.0;
    //       dsraddz = 0.0d0
    dsraddz = 0.0;
    //
    //
    //! ion section:
    //        xni     = avo * ytot1 * den
    xni = ((avo * ytot1) * den);
    //        dxnidd  = avo * ytot1
    dxnidd = (avo * ytot1);
    //        dxnida  = -xni * ytot1
    dxnida = ((-xni) * ytot1);
    //
    //        pion    = xni * kt
    pion = (xni * kt);
    //        dpiondd = dxnidd * kt
    dpiondd = (dxnidd * kt);
    //        dpiondt = xni * kerg
    dpiondt = (xni * kerg);
    //        dpionda = dxnida * kt
    dpionda = (dxnida * kt);
    //        dpiondz = 0.0d0
    dpiondz = 0.0;
    //
    //        eion    = 1.5d0 * pion*deni
    eion = ((1.5 * pion) * deni);
    //        deiondd = (1.5d0 * dpiondd - eion)*deni
    deiondd = (((1.5 * dpiondd) - eion) * deni);
    //        deiondt = 1.5d0 * dpiondt*deni
    deiondt = ((1.5 * dpiondt) * deni);
    //        deionda = 1.5d0 * dpionda*deni
    deionda = ((1.5 * dpionda) * deni);
    //        deiondz = 0.0d0
    deiondz = 0.0;
    //
    //
    //! sackur-tetrode equation for the ion entropy of
    //! a single ideal gas characterized by abar
    //        x       = abar*abar*sqrt(abar) * deni/avo
    x = ((((abar * abar) * std::sqrt(abar)) * deni) / avo);
    //        s       = sioncon * temp
    s = (sioncon * temp);
    //        z       = x * s * sqrt(s)
    z = ((x * s) * std::sqrt(s));
    //        y       = log(z)
    y = std::log(z);
    //
    //!        y       = 1.0d0/(abar*kt)
    //!        yy      = y * sqrt(y)
    //!        z       = xni * sifac * yy
    //!        etaion  = log(z)
    //
    //
    //        sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
    sion = ((((pion * deni) + eion) * tempi) + ((kergavo * ytot1) * y));
    //        dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
    //                   - kergavo * deni * ytot1
    dsiondd =
        (((((dpiondd * deni) - ((pion * deni) * deni)) + deiondd) * tempi) - ((kergavo * deni) * ytot1));
    //        dsiondt = (dpiondt*deni + deiondt)*tempi - &
    //                  (pion*deni + eion) * tempi*tempi &
    //                  + 1.5d0 * kergavo * tempi*ytot1
    dsiondt = (((((dpiondt * deni) + deiondt) * tempi) - ((((pion * deni) + eion) * tempi) * tempi)) +
               (((1.5 * kergavo) * tempi) * ytot1));
    //        x       = avo*kerg/abar
    x = ((avo * kerg) / abar);
    //        dsionda = (dpionda*deni + deionda)*tempi &
    //                  + kergavo*ytot1*ytot1* (2.5d0 - y)
    dsionda = ((((dpionda * deni) + deionda) * tempi) + (((kergavo * ytot1) * ytot1) * (2.5 - y)));
    //        dsiondz = 0.0d0
    dsiondz = 0.0;
    //
    //
    //
    //! electron-positron section:
    //
    //
    //! assume complete ionization
    //        xnem    = xni * zbar
    xnem = (xni * zbar);
    //
    //
    //! enter the table with ye*den
    //        din = ye*den
    din = (ye * den);
    //
    //
    //! bomb proof the input
    //        if (temp .gt. t(jmax)) then
    if ((temp > table.t[jmax])) {
        //         write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
        //         write(6,*) 'temp too hot, off grid'
        //         write(6,*) 'setting eosfail to true and returning'
        //         eosfail = .true.
        throw std::out_of_range("Helmholtz state is outside the table");
        //         return
        //        end if
    }
    //        if (temp .lt. t(1)) then
    if ((temp < table.t[1])) {
        //         write(6,01) 'temp=',temp,' t(1)=',t(1)
        //         write(6,*) 'temp too cold, off grid'
        //         write(6,*) 'setting eosfail to true and returning'
        //         eosfail = .true.
        throw std::out_of_range("Helmholtz state is outside the table");
        //         return
        //        end if
    }
    //        if (din  .gt. d(imax)) then
    if ((din > table.d[imax])) {
        //         write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
        //         write(6,*) 'ye*den too big, off grid'
        //         write(6,*) 'setting eosfail to true and returning'
        //         eosfail = .true.
        throw std::out_of_range("Helmholtz state is outside the table");
        //         return
        //        end if
    }
    //        if (din  .lt. d(1)) then
    if ((din < table.d[1])) {
        //         write(6,01) 'ye*den=',din,' d(1)=',d(1)
        //         write(6,*) 'ye*den too small, off grid'
        //         write(6,*) 'setting eosfail to true and returning'
        //         eosfail = .true.
        throw std::out_of_range("Helmholtz state is outside the table");
        //         return
        //        end if
    }
    //
    //! hash locate this temperature and density
    //        jat = int((log10(temp) - tlo)*tstpi) + 1
    jat = (static_cast<int>(((std::log10(temp) - tlo) * tstpi)) + 1);
    //        jat = max(1,min(jat,jmax-1))
    jat = std::max(1, std::min(jat, (jmax - 1)));
    //        iat = int((log10(din) - dlo)*dstpi) + 1
    iat = (static_cast<int>(((std::log10(din) - dlo) * dstpi)) + 1);
    //        iat = max(1,min(iat,imax-1))
    iat = std::max(1, std::min(iat, (imax - 1)));
    //
    //
    //! access the table locations only once
    //        fi(1)  = f(iat,jat)
    fi[1] = table.f(iat, jat);
    //        fi(2)  = f(iat+1,jat)
    fi[2] = table.f((iat + 1), jat);
    //        fi(3)  = f(iat,jat+1)
    fi[3] = table.f(iat, (jat + 1));
    //        fi(4)  = f(iat+1,jat+1)
    fi[4] = table.f((iat + 1), (jat + 1));
    //        fi(5)  = ft(iat,jat)
    fi[5] = table.ft(iat, jat);
    //        fi(6)  = ft(iat+1,jat)
    fi[6] = table.ft((iat + 1), jat);
    //        fi(7)  = ft(iat,jat+1)
    fi[7] = table.ft(iat, (jat + 1));
    //        fi(8)  = ft(iat+1,jat+1)
    fi[8] = table.ft((iat + 1), (jat + 1));
    //        fi(9)  = ftt(iat,jat)
    fi[9] = table.ftt(iat, jat);
    //        fi(10) = ftt(iat+1,jat)
    fi[10] = table.ftt((iat + 1), jat);
    //        fi(11) = ftt(iat,jat+1)
    fi[11] = table.ftt(iat, (jat + 1));
    //        fi(12) = ftt(iat+1,jat+1)
    fi[12] = table.ftt((iat + 1), (jat + 1));
    //        fi(13) = fd(iat,jat)
    fi[13] = table.fd(iat, jat);
    //        fi(14) = fd(iat+1,jat)
    fi[14] = table.fd((iat + 1), jat);
    //        fi(15) = fd(iat,jat+1)
    fi[15] = table.fd(iat, (jat + 1));
    //        fi(16) = fd(iat+1,jat+1)
    fi[16] = table.fd((iat + 1), (jat + 1));
    //        fi(17) = fdd(iat,jat)
    fi[17] = table.fdd(iat, jat);
    //        fi(18) = fdd(iat+1,jat)
    fi[18] = table.fdd((iat + 1), jat);
    //        fi(19) = fdd(iat,jat+1)
    fi[19] = table.fdd(iat, (jat + 1));
    //        fi(20) = fdd(iat+1,jat+1)
    fi[20] = table.fdd((iat + 1), (jat + 1));
    //        fi(21) = fdt(iat,jat)
    fi[21] = table.fdt(iat, jat);
    //        fi(22) = fdt(iat+1,jat)
    fi[22] = table.fdt((iat + 1), jat);
    //        fi(23) = fdt(iat,jat+1)
    fi[23] = table.fdt(iat, (jat + 1));
    //        fi(24) = fdt(iat+1,jat+1)
    fi[24] = table.fdt((iat + 1), (jat + 1));
    //        fi(25) = fddt(iat,jat)
    fi[25] = table.fddt(iat, jat);
    //        fi(26) = fddt(iat+1,jat)
    fi[26] = table.fddt((iat + 1), jat);
    //        fi(27) = fddt(iat,jat+1)
    fi[27] = table.fddt(iat, (jat + 1));
    //        fi(28) = fddt(iat+1,jat+1)
    fi[28] = table.fddt((iat + 1), (jat + 1));
    //        fi(29) = fdtt(iat,jat)
    fi[29] = table.fdtt(iat, jat);
    //        fi(30) = fdtt(iat+1,jat)
    fi[30] = table.fdtt((iat + 1), jat);
    //        fi(31) = fdtt(iat,jat+1)
    fi[31] = table.fdtt(iat, (jat + 1));
    //        fi(32) = fdtt(iat+1,jat+1)
    fi[32] = table.fdtt((iat + 1), (jat + 1));
    //        fi(33) = fddtt(iat,jat)
    fi[33] = table.fddtt(iat, jat);
    //        fi(34) = fddtt(iat+1,jat)
    fi[34] = table.fddtt((iat + 1), jat);
    //        fi(35) = fddtt(iat,jat+1)
    fi[35] = table.fddtt(iat, (jat + 1));
    //        fi(36) = fddtt(iat+1,jat+1)
    fi[36] = table.fddtt((iat + 1), (jat + 1));
    //
    //
    //! various differences
    //        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
    xt = std::max(((temp - table.t[jat]) * table.dti_sav[jat]), 0.0);
    //        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
    xd = std::max(((din - table.d[iat]) * table.ddi_sav[iat]), 0.0);
    //        mxt = 1.0d0 - xt
    mxt = (1.0 - xt);
    //        mxd = 1.0d0 - xd
    mxd = (1.0 - xd);
    //
    //! the six density and six temperature basis functions
    //        si0t =   psi0(xt)
    si0t = psi0(xt);
    //        si1t =   psi1(xt)*dt_sav(jat)
    si1t = (psi1(xt) * table.dt_sav[jat]);
    //        si2t =   psi2(xt)*dt2_sav(jat)
    si2t = (psi2(xt) * table.dt2_sav[jat]);
    //
    //        si0mt =  psi0(mxt)
    si0mt = psi0(mxt);
    //        si1mt = -psi1(mxt)*dt_sav(jat)
    si1mt = ((-psi1(mxt)) * table.dt_sav[jat]);
    //        si2mt =  psi2(mxt)*dt2_sav(jat)
    si2mt = (psi2(mxt) * table.dt2_sav[jat]);
    //
    //        si0d =   psi0(xd)
    si0d = psi0(xd);
    //        si1d =   psi1(xd)*dd_sav(iat)
    si1d = (psi1(xd) * table.dd_sav[iat]);
    //        si2d =   psi2(xd)*dd2_sav(iat)
    si2d = (psi2(xd) * table.dd2_sav[iat]);
    //
    //        si0md =  psi0(mxd)
    si0md = psi0(mxd);
    //        si1md = -psi1(mxd)*dd_sav(iat)
    si1md = ((-psi1(mxd)) * table.dd_sav[iat]);
    //        si2md =  psi2(mxd)*dd2_sav(iat)
    si2md = (psi2(mxd) * table.dd2_sav[iat]);
    //
    //! derivatives of the weight functions
    //        dsi0t =   dpsi0(xt)*dti_sav(jat)
    dsi0t = (dpsi0(xt) * table.dti_sav[jat]);
    //        dsi1t =   dpsi1(xt)
    dsi1t = dpsi1(xt);
    //        dsi2t =   dpsi2(xt)*dt_sav(jat)
    dsi2t = (dpsi2(xt) * table.dt_sav[jat]);
    //
    //        dsi0mt = -dpsi0(mxt)*dti_sav(jat)
    dsi0mt = ((-dpsi0(mxt)) * table.dti_sav[jat]);
    //        dsi1mt =  dpsi1(mxt)
    dsi1mt = dpsi1(mxt);
    //        dsi2mt = -dpsi2(mxt)*dt_sav(jat)
    dsi2mt = ((-dpsi2(mxt)) * table.dt_sav[jat]);
    //
    //        dsi0d =   dpsi0(xd)*ddi_sav(iat)
    dsi0d = (dpsi0(xd) * table.ddi_sav[iat]);
    //        dsi1d =   dpsi1(xd)
    dsi1d = dpsi1(xd);
    //        dsi2d =   dpsi2(xd)*dd_sav(iat)
    dsi2d = (dpsi2(xd) * table.dd_sav[iat]);
    //
    //        dsi0md = -dpsi0(mxd)*ddi_sav(iat)
    dsi0md = ((-dpsi0(mxd)) * table.ddi_sav[iat]);
    //        dsi1md =  dpsi1(mxd)
    dsi1md = dpsi1(mxd);
    //        dsi2md = -dpsi2(mxd)*dd_sav(iat)
    dsi2md = ((-dpsi2(mxd)) * table.dd_sav[iat]);
    //
    //! second derivatives of the weight functions
    //        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
    ddsi0t = (ddpsi0(xt) * table.dt2i_sav[jat]);
    //        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
    ddsi1t = (ddpsi1(xt) * table.dti_sav[jat]);
    //        ddsi2t =   ddpsi2(xt)
    ddsi2t = ddpsi2(xt);
    //
    //        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
    ddsi0mt = (ddpsi0(mxt) * table.dt2i_sav[jat]);
    //        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
    ddsi1mt = ((-ddpsi1(mxt)) * table.dti_sav[jat]);
    //        ddsi2mt =  ddpsi2(mxt)
    ddsi2mt = ddpsi2(mxt);
    //
    //!        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
    //!        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
    //!        ddsi2d =   ddpsi2(xd)
    //
    //!        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
    //!        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
    //!        ddsi2md =  ddpsi2(mxd)
    //
    //
    //! the free energy
    //        free  = h5(iat,jat, &
    //                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
    //                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
    free = h5(si0t, si1t, si2t, si0mt, si1mt, si2mt, si0d, si1d, si2d, si0md, si1md, si2md);
    //
    //! derivative with respect to density
    //        df_d  = h5(iat,jat, &
    //                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
    //                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)
    df_d = h5(si0t, si1t, si2t, si0mt, si1mt, si2mt, dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md);
    //
    //
    //! derivative with respect to temperature
    //        df_t = h5(iat,jat, &
    //                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
    //                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
    df_t = h5(dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, si0d, si1d, si2d, si0md, si1md, si2md);
    //
    //! derivative with respect to density**2
    //!        df_dd = h5(iat,jat,
    //!     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
    //!     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)
    //
    //! derivative with respect to temperature**2
    //        df_tt = h5(iat,jat, &
    //              ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
    //                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
    df_tt = h5(ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, si0d, si1d, si2d, si0md, si1md, si2md);
    //
    //! derivative with respect to temperature and density
    //        df_dt = h5(iat,jat, &
    //                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
    //                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)
    df_dt = h5(dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md);
    //
    //
    //
    //! now get the pressure derivative with density, chemical potential, and
    //! electron positron number densities
    //! get the interpolation weight functions
    //        si0t   =  xpsi0(xt)
    si0t = xpsi0(xt);
    //        si1t   =  xpsi1(xt)*dt_sav(jat)
    si1t = (xpsi1(xt) * table.dt_sav[jat]);
    //
    //        si0mt  =  xpsi0(mxt)
    si0mt = xpsi0(mxt);
    //        si1mt  =  -xpsi1(mxt)*dt_sav(jat)
    si1mt = ((-xpsi1(mxt)) * table.dt_sav[jat]);
    //
    //        si0d   =  xpsi0(xd)
    si0d = xpsi0(xd);
    //        si1d   =  xpsi1(xd)*dd_sav(iat)
    si1d = (xpsi1(xd) * table.dd_sav[iat]);
    //
    //        si0md  =  xpsi0(mxd)
    si0md = xpsi0(mxd);
    //        si1md  =  -xpsi1(mxd)*dd_sav(iat)
    si1md = ((-xpsi1(mxd)) * table.dd_sav[iat]);
    //
    //
    //! derivatives of weight functions
    //        dsi0t  = xdpsi0(xt)*dti_sav(jat)
    dsi0t = (xdpsi0(xt) * table.dti_sav[jat]);
    //        dsi1t  = xdpsi1(xt)
    dsi1t = xdpsi1(xt);
    //
    //        dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
    dsi0mt = ((-xdpsi0(mxt)) * table.dti_sav[jat]);
    //        dsi1mt = xdpsi1(mxt)
    dsi1mt = xdpsi1(mxt);
    //
    //        dsi0d  = xdpsi0(xd)*ddi_sav(iat)
    dsi0d = (xdpsi0(xd) * table.ddi_sav[iat]);
    //        dsi1d  = xdpsi1(xd)
    dsi1d = xdpsi1(xd);
    //
    //        dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
    dsi0md = ((-xdpsi0(mxd)) * table.ddi_sav[iat]);
    //        dsi1md = xdpsi1(mxd)
    dsi1md = xdpsi1(mxd);
    //
    //
    //! look in the pressure derivative only once
    //        fi(1)  = dpdf(iat,jat)
    fi[1] = table.dpdf(iat, jat);
    //        fi(2)  = dpdf(iat+1,jat)
    fi[2] = table.dpdf((iat + 1), jat);
    //        fi(3)  = dpdf(iat,jat+1)
    fi[3] = table.dpdf(iat, (jat + 1));
    //        fi(4)  = dpdf(iat+1,jat+1)
    fi[4] = table.dpdf((iat + 1), (jat + 1));
    //        fi(5)  = dpdft(iat,jat)
    fi[5] = table.dpdft(iat, jat);
    //        fi(6)  = dpdft(iat+1,jat)
    fi[6] = table.dpdft((iat + 1), jat);
    //        fi(7)  = dpdft(iat,jat+1)
    fi[7] = table.dpdft(iat, (jat + 1));
    //        fi(8)  = dpdft(iat+1,jat+1)
    fi[8] = table.dpdft((iat + 1), (jat + 1));
    //        fi(9)  = dpdfd(iat,jat)
    fi[9] = table.dpdfd(iat, jat);
    //        fi(10) = dpdfd(iat+1,jat)
    fi[10] = table.dpdfd((iat + 1), jat);
    //        fi(11) = dpdfd(iat,jat+1)
    fi[11] = table.dpdfd(iat, (jat + 1));
    //        fi(12) = dpdfd(iat+1,jat+1)
    fi[12] = table.dpdfd((iat + 1), (jat + 1));
    //        fi(13) = dpdfdt(iat,jat)
    fi[13] = table.dpdfdt(iat, jat);
    //        fi(14) = dpdfdt(iat+1,jat)
    fi[14] = table.dpdfdt((iat + 1), jat);
    //        fi(15) = dpdfdt(iat,jat+1)
    fi[15] = table.dpdfdt(iat, (jat + 1));
    //        fi(16) = dpdfdt(iat+1,jat+1)
    fi[16] = table.dpdfdt((iat + 1), (jat + 1));
    //
    //! pressure derivative with density
    //        dpepdd  = h3(iat,jat, &
    //                       si0t,   si1t,   si0mt,   si1mt, &
    //                       si0d,   si1d,   si0md,   si1md)
    dpepdd = h3(si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md);
    //        dpepdd  = max(ye * dpepdd,1.0d-30)
    dpepdd = std::max((ye * dpepdd), 1e-30);
    //
    //
    //
    //! look in the electron chemical potential table only once
    //        fi(1)  = ef(iat,jat)
    fi[1] = table.ef(iat, jat);
    //        fi(2)  = ef(iat+1,jat)
    fi[2] = table.ef((iat + 1), jat);
    //        fi(3)  = ef(iat,jat+1)
    fi[3] = table.ef(iat, (jat + 1));
    //        fi(4)  = ef(iat+1,jat+1)
    fi[4] = table.ef((iat + 1), (jat + 1));
    //        fi(5)  = eft(iat,jat)
    fi[5] = table.eft(iat, jat);
    //        fi(6)  = eft(iat+1,jat)
    fi[6] = table.eft((iat + 1), jat);
    //        fi(7)  = eft(iat,jat+1)
    fi[7] = table.eft(iat, (jat + 1));
    //        fi(8)  = eft(iat+1,jat+1)
    fi[8] = table.eft((iat + 1), (jat + 1));
    //        fi(9)  = efd(iat,jat)
    fi[9] = table.efd(iat, jat);
    //        fi(10) = efd(iat+1,jat)
    fi[10] = table.efd((iat + 1), jat);
    //        fi(11) = efd(iat,jat+1)
    fi[11] = table.efd(iat, (jat + 1));
    //        fi(12) = efd(iat+1,jat+1)
    fi[12] = table.efd((iat + 1), (jat + 1));
    //        fi(13) = efdt(iat,jat)
    fi[13] = table.efdt(iat, jat);
    //        fi(14) = efdt(iat+1,jat)
    fi[14] = table.efdt((iat + 1), jat);
    //        fi(15) = efdt(iat,jat+1)
    fi[15] = table.efdt(iat, (jat + 1));
    //        fi(16) = efdt(iat+1,jat+1)
    fi[16] = table.efdt((iat + 1), (jat + 1));
    //
    //
    //! electron chemical potential etaele
    //        etaele  = h3(iat,jat, &
    //                     si0t,   si1t,   si0mt,   si1mt, &
    //                     si0d,   si1d,   si0md,   si1md)
    etaele = h3(si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md);
    //
    //
    //! derivative with respect to density
    //        x       = h3(iat,jat, &
    //                     si0t,   si1t,   si0mt,   si1mt, &
    //                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
    x = h3(si0t, si1t, si0mt, si1mt, dsi0d, dsi1d, dsi0md, dsi1md);
    //        detadd  = ye * x
    detadd = (ye * x);
    //
    //! derivative with respect to temperature
    //        detadt  = h3(iat,jat, &
    //                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
    //                     si0d,   si1d,   si0md,   si1md)
    detadt = h3(dsi0t, dsi1t, dsi0mt, dsi1mt, si0d, si1d, si0md, si1md);
    //
    //! derivative with respect to abar and zbar
    //       detada = -x * din * ytot1
    detada = (((-x) * din) * ytot1);
    //       detadz =  x * den * ytot1
    detadz = ((x * den) * ytot1);
    //
    //
    //
    //! look in the number density table only once
    //        fi(1)  = xf(iat,jat)
    fi[1] = table.xf(iat, jat);
    //        fi(2)  = xf(iat+1,jat)
    fi[2] = table.xf((iat + 1), jat);
    //        fi(3)  = xf(iat,jat+1)
    fi[3] = table.xf(iat, (jat + 1));
    //        fi(4)  = xf(iat+1,jat+1)
    fi[4] = table.xf((iat + 1), (jat + 1));
    //        fi(5)  = xft(iat,jat)
    fi[5] = table.xft(iat, jat);
    //        fi(6)  = xft(iat+1,jat)
    fi[6] = table.xft((iat + 1), jat);
    //        fi(7)  = xft(iat,jat+1)
    fi[7] = table.xft(iat, (jat + 1));
    //        fi(8)  = xft(iat+1,jat+1)
    fi[8] = table.xft((iat + 1), (jat + 1));
    //        fi(9)  = xfd(iat,jat)
    fi[9] = table.xfd(iat, jat);
    //        fi(10) = xfd(iat+1,jat)
    fi[10] = table.xfd((iat + 1), jat);
    //        fi(11) = xfd(iat,jat+1)
    fi[11] = table.xfd(iat, (jat + 1));
    //        fi(12) = xfd(iat+1,jat+1)
    fi[12] = table.xfd((iat + 1), (jat + 1));
    //        fi(13) = xfdt(iat,jat)
    fi[13] = table.xfdt(iat, jat);
    //        fi(14) = xfdt(iat+1,jat)
    fi[14] = table.xfdt((iat + 1), jat);
    //        fi(15) = xfdt(iat,jat+1)
    fi[15] = table.xfdt(iat, (jat + 1));
    //        fi(16) = xfdt(iat+1,jat+1)
    fi[16] = table.xfdt((iat + 1), (jat + 1));
    //
    //! electron + positron number densities
    //       xnefer   = h3(iat,jat, &
    //                     si0t,   si1t,   si0mt,   si1mt, &
    //                     si0d,   si1d,   si0md,   si1md)
    xnefer = h3(si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md);
    //
    //! derivative with respect to density
    //       x        = h3(iat,jat, &
    //                     si0t,   si1t,   si0mt,   si1mt, &
    //                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
    x = h3(si0t, si1t, si0mt, si1mt, dsi0d, dsi1d, dsi0md, dsi1md);
    //       x = max(x,1.0d-30)
    x = std::max(x, 1e-30);
    //       dxnedd   = ye * x
    dxnedd = (ye * x);
    //
    //! derivative with respect to temperature
    //       dxnedt   = h3(iat,jat, &
    //                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
    //                     si0d,   si1d,   si0md,   si1md)
    dxnedt = h3(dsi0t, dsi1t, dsi0mt, dsi1mt, si0d, si1d, si0md, si1md);
    //
    //! derivative with respect to abar and zbar
    //       dxneda = -x * din * ytot1
    dxneda = (((-x) * din) * ytot1);
    //       dxnedz =  x  * den * ytot1
    dxnedz = ((x * den) * ytot1);
    //
    //
    //! the desired electron-positron thermodynamic quantities
    //
    //! dpepdd at high temperatures and low densities is below the
    //! floating point limit of the subtraction of two large terms.
    //! since dpresdd doesn't enter the maxwell relations at all, use the
    //! bicubic interpolation done above instead of the formally correct expression
    //        x       = din * din
    x = (din * din);
    //        pele    = x * df_d
    pele = (x * df_d);
    //        dpepdt  = x * df_dt
    dpepdt = (x * df_dt);
    //!        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
    //        s       = dpepdd/ye - 2.0d0 * din * df_d
    s = ((dpepdd / ye) - ((2.0 * din) * df_d));
    //        dpepda  = -ytot1 * (2.0d0 * pele + s * din)
    dpepda = ((-ytot1) * ((2.0 * pele) + (s * din)));
    //        dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)
    dpepdz = ((den * ytot1) * (((2.0 * din) * df_d) + s));
    //
    //
    //        x       = ye * ye
    x = (ye * ye);
    //        sele    = -df_t * ye
    sele = ((-df_t) * ye);
    //        dsepdt  = -df_tt * ye
    dsepdt = ((-df_tt) * ye);
    //        dsepdd  = -df_dt * x
    dsepdd = ((-df_dt) * x);
    //        dsepda  = ytot1 * (ye * df_dt * din - sele)
    dsepda = (ytot1 * (((ye * df_dt) * din) - sele));
    //        dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)
    dsepdz = ((-ytot1) * (((ye * df_dt) * den) + df_t));
    //
    //
    //        eele    = ye*free + temp * sele
    eele = ((ye * free) + (temp * sele));
    //        deepdt  = temp * dsepdt
    deepdt = (temp * dsepdt);
    //        deepdd  = x * df_d + temp * dsepdd
    deepdd = ((x * df_d) + (temp * dsepdd));
    //        deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
    deepda = ((((-ye) * ytot1) * (free + (df_d * din))) + (temp * dsepda));
    //        deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz
    deepdz = ((ytot1 * (free + ((ye * df_d) * den))) + (temp * dsepdz));
    //
    //
    //
    //
    //! coulomb section:
    //
    //! uniform background corrections only
    //! from yakovlev & shalybkov 1989
    //! lami is the average ion seperation
    //! plasg is the plasma coupling parameter
    //
    //        z        = forth * pi
    z = (forth * pi);
    //        s        = z * xni
    s = (z * xni);
    //        dsdd     = z * dxnidd
    dsdd = (z * dxnidd);
    //        dsda     = z * dxnida
    dsda = (z * dxnida);
    //
    //        lami     = 1.0d0/s**third
    lami = (1.0 / std::pow(s, third));
    //        inv_lami = 1.0d0/lami
    inv_lami = (1.0 / lami);
    //        z        = -third * lami
    z = ((-third) * lami);
    //        lamidd   = z * dsdd/s
    lamidd = ((z * dsdd) / s);
    //        lamida   = z * dsda/s
    lamida = ((z * dsda) / s);
    //
    //        plasg    = zbar*zbar*esqu*ktinv*inv_lami
    plasg = ((((zbar * zbar) * esqu) * ktinv) * inv_lami);
    //        z        = -plasg * inv_lami
    z = ((-plasg) * inv_lami);
    //        plasgdd  = z * lamidd
    plasgdd = (z * lamidd);
    //        plasgda  = z * lamida
    plasgda = (z * lamida);
    //        plasgdt  = -plasg*ktinv * kerg
    plasgdt = (((-plasg) * ktinv) * kerg);
    //        plasgdz  = 2.0d0 * plasg/zbar
    plasgdz = ((2.0 * plasg) / zbar);
    //
    //
    //! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    //        if (plasg .ge. 1.0) then
    if ((plasg >= 1.0)) {
        //         x        = plasg**(0.25d0)
        x = std::pow(plasg, 0.25);
        //         y        = avo * ytot1 * kerg
        y = ((avo * ytot1) * kerg);
        //         ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
        ecoul = ((y * temp) * ((((a1 * plasg) + (b1 * x)) + (c1 / x)) + d1));
        //         pcoul    = third * den * ecoul
        pcoul = ((third * den) * ecoul);
        //         scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
        //                    + d1 * (log(plasg) - 1.0d0) - e1)
        scoul = ((-y) * (((((3.0 * b1) * x) - ((5.0 * c1) / x)) + (d1 * (std::log(plasg) - 1.0))) - e1));
        //
        //         y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
        y = (((avo * ytot1) * kt) * (a1 + ((0.25 / plasg) * ((b1 * x) - (c1 / x)))));
        //         decouldd = y * plasgdd
        decouldd = (y * plasgdd);
        //         decouldt = y * plasgdt + ecoul/temp
        decouldt = ((y * plasgdt) + (ecoul / temp));
        //         decoulda = y * plasgda - ecoul/abar
        decoulda = ((y * plasgda) - (ecoul / abar));
        //         decouldz = y * plasgdz
        decouldz = (y * plasgdz);
        //
        //         y        = third * den
        y = (third * den);
        //         dpcouldd = third * ecoul + y*decouldd
        dpcouldd = ((third * ecoul) + (y * decouldd));
        //         dpcouldt = y * decouldt
        dpcouldt = (y * decouldt);
        //         dpcoulda = y * decoulda
        dpcoulda = (y * decoulda);
        //         dpcouldz = y * decouldz
        dpcouldz = (y * decouldz);
        //
        //
        //         y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
        y = ((((-avo) * kerg) / (abar * plasg)) * ((((0.75 * b1) * x) + ((1.25 * c1) / x)) + d1));
        //         dscouldd = y * plasgdd
        dscouldd = (y * plasgdd);
        //         dscouldt = y * plasgdt
        dscouldt = (y * plasgdt);
        //         dscoulda = y * plasgda - scoul/abar
        dscoulda = ((y * plasgda) - (scoul / abar));
        //         dscouldz = y * plasgdz
        dscouldz = (y * plasgdz);
        //
        //
        //! yakovlev & shalybkov 1989 equations 102, 103, 104
        //        else if (plasg .lt. 1.0) then
    } else if ((plasg < 1.0)) {
        //         x        = plasg*sqrt(plasg)
        x = (plasg * std::sqrt(plasg));
        //         y        = plasg**b2
        y = std::pow(plasg, b2);
        //         z        = c2 * x - third * a2 * y
        z = ((c2 * x) - ((third * a2) * y));
        //         pcoul    = -pion * z
        pcoul = ((-pion) * z);
        //         ecoul    = 3.0d0 * pcoul/den
        ecoul = ((3.0 * pcoul) / den);
        //         scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)
        scoul = ((((-avo) / abar) * kerg) * ((c2 * x) - (((a2 * (b2 - 1.0)) / b2) * y)));
        //
        //         s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
        s = ((((1.5 * c2) * x) / plasg) - ((((third * a2) * b2) * y) / plasg));
        //         dpcouldd = -dpiondd*z - pion*s*plasgdd
        dpcouldd = (((-dpiondd) * z) - ((pion * s) * plasgdd));
        //         dpcouldt = -dpiondt*z - pion*s*plasgdt
        dpcouldt = (((-dpiondt) * z) - ((pion * s) * plasgdt));
        //         dpcoulda = -dpionda*z - pion*s*plasgda
        dpcoulda = (((-dpionda) * z) - ((pion * s) * plasgda));
        //         dpcouldz = -dpiondz*z - pion*s*plasgdz
        dpcouldz = (((-dpiondz) * z) - ((pion * s) * plasgdz));
        //
        //         s        = 3.0d0/den
        s = (3.0 / den);
        //         decouldd = s * dpcouldd - ecoul/den
        decouldd = ((s * dpcouldd) - (ecoul / den));
        //         decouldt = s * dpcouldt
        decouldt = (s * dpcouldt);
        //         decoulda = s * dpcoulda
        decoulda = (s * dpcoulda);
        //         decouldz = s * dpcouldz
        decouldz = (s * dpcouldz);
        //
        //         s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
        s = ((((-avo) * kerg) / (abar * plasg)) * (((1.5 * c2) * x) - ((a2 * (b2 - 1.0)) * y)));
        //         dscouldd = s * plasgdd
        dscouldd = (s * plasgdd);
        //         dscouldt = s * plasgdt
        dscouldt = (s * plasgdt);
        //         dscoulda = s * plasgda - scoul/abar
        dscoulda = ((s * plasgda) - (scoul / abar));
        //         dscouldz = s * plasgdz
        dscouldz = (s * plasgdz);
        //        end if
    }
    //
    //
    //! bomb proof
    //        x   = prad + pion + pele + pcoul
    x = (((prad + pion) + pele) + pcoul);
    //        y   = erad + eion + eele + ecoul
    y = (((erad + eion) + eele) + ecoul);
    //        z   = srad + sion + sele + scoul
    z = (((srad + sion) + sele) + scoul);
    //
    //!        write(6,*) x,y,z
    //!        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
    //        if (x .le. 0.0 .or. y .le. 0.0) then
    if (((x <= 0.0) || (y <= 0.0))) {
        //!        if (x .le. 0.0) then
        //
        //!         write(6,*)
        //!         write(6,*) 'coulomb corrections are causing a negative pressure'
        //!         write(6,*) 'setting all coulomb corrections to zero'
        //!         write(6,*)
        //
        //         pcoul    = 0.0d0
        pcoul = 0.0;
        //         dpcouldd = 0.0d0
        dpcouldd = 0.0;
        //         dpcouldt = 0.0d0
        dpcouldt = 0.0;
        //         dpcoulda = 0.0d0
        dpcoulda = 0.0;
        //         dpcouldz = 0.0d0
        dpcouldz = 0.0;
        //         ecoul    = 0.0d0
        ecoul = 0.0;
        //         decouldd = 0.0d0
        decouldd = 0.0;
        //         decouldt = 0.0d0
        decouldt = 0.0;
        //         decoulda = 0.0d0
        decoulda = 0.0;
        //         decouldz = 0.0d0
        decouldz = 0.0;
        //         scoul    = 0.0d0
        scoul = 0.0;
        //         dscouldd = 0.0d0
        dscouldd = 0.0;
        //         dscouldt = 0.0d0
        dscouldt = 0.0;
        //         dscoulda = 0.0d0
        dscoulda = 0.0;
        //         dscouldz = 0.0d0
        dscouldz = 0.0;
        //        end if
    }
    //
    //
    //! sum all the gas components
    //       pgas    = pion + pele + pcoul
    pgas = ((pion + pele) + pcoul);
    //       egas    = eion + eele + ecoul
    egas = ((eion + eele) + ecoul);
    //       sgas    = sion + sele + scoul
    sgas = ((sion + sele) + scoul);
    //
    //       dpgasdd = dpiondd + dpepdd + dpcouldd
    dpgasdd = ((dpiondd + dpepdd) + dpcouldd);
    //       dpgasdt = dpiondt + dpepdt + dpcouldt
    dpgasdt = ((dpiondt + dpepdt) + dpcouldt);
    //       dpgasda = dpionda + dpepda + dpcoulda
    dpgasda = ((dpionda + dpepda) + dpcoulda);
    //       dpgasdz = dpiondz + dpepdz + dpcouldz
    dpgasdz = ((dpiondz + dpepdz) + dpcouldz);
    //
    //       degasdd = deiondd + deepdd + decouldd
    degasdd = ((deiondd + deepdd) + decouldd);
    //       degasdt = deiondt + deepdt + decouldt
    degasdt = ((deiondt + deepdt) + decouldt);
    //       degasda = deionda + deepda + decoulda
    degasda = ((deionda + deepda) + decoulda);
    //       degasdz = deiondz + deepdz + decouldz
    degasdz = ((deiondz + deepdz) + decouldz);
    //
    //       dsgasdd = dsiondd + dsepdd + dscouldd
    dsgasdd = ((dsiondd + dsepdd) + dscouldd);
    //       dsgasdt = dsiondt + dsepdt + dscouldt
    dsgasdt = ((dsiondt + dsepdt) + dscouldt);
    //       dsgasda = dsionda + dsepda + dscoulda
    dsgasda = ((dsionda + dsepda) + dscoulda);
    //       dsgasdz = dsiondz + dsepdz + dscouldz
    dsgasdz = ((dsiondz + dsepdz) + dscouldz);
    //
    //
    //
    //
    //! add in radiation to get the total
    //       pres    = prad + pgas
    pres = (prad + pgas);
    //       ener    = erad + egas
    ener = (erad + egas);
    //       entr    = srad + sgas
    entr = (srad + sgas);
    //
    //       dpresdd = dpraddd + dpgasdd
    dpresdd = (dpraddd + dpgasdd);
    //       dpresdt = dpraddt + dpgasdt
    dpresdt = (dpraddt + dpgasdt);
    //       dpresda = dpradda + dpgasda
    dpresda = (dpradda + dpgasda);
    //       dpresdz = dpraddz + dpgasdz
    dpresdz = (dpraddz + dpgasdz);
    //
    //       denerdd = deraddd + degasdd
    denerdd = (deraddd + degasdd);
    //       denerdt = deraddt + degasdt
    denerdt = (deraddt + degasdt);
    //       denerda = deradda + degasda
    denerda = (deradda + degasda);
    //       denerdz = deraddz + degasdz
    denerdz = (deraddz + degasdz);
    //
    //       dentrdd = dsraddd + dsgasdd
    dentrdd = (dsraddd + dsgasdd);
    //       dentrdt = dsraddt + dsgasdt
    dentrdt = (dsraddt + dsgasdt);
    //       dentrda = dsradda + dsgasda
    dentrda = (dsradda + dsgasda);
    //       dentrdz = dsraddz + dsgasdz
    dentrdz = (dsraddz + dsgasdz);
    //
    //
    //! for the gas
    //! the temperature and density exponents (c&g 9.81 9.82)
    //! the specific heat at constant volume (c&g 9.92)
    //! the third adiabatic exponent (c&g 9.93)
    //! the first adiabatic exponent (c&g 9.97)
    //! the second adiabatic exponent (c&g 9.105)
    //! the specific heat at constant pressure (c&g 9.98)
    //! and relativistic formula for the sound speed (c&g 14.29)
    //
    //       zz        = pgas*deni
    zz = (pgas * deni);
    //       zzi       = den/pgas
    zzi = (den / pgas);
    //       chit_gas  = temp/pgas * dpgasdt
    chit_gas = ((temp / pgas) * dpgasdt);
    //       chid_gas  = dpgasdd*zzi
    chid_gas = (dpgasdd * zzi);
    //       cv_gas    = degasdt
    cv_gas = degasdt;
    //       x         = zz * chit_gas/(temp * cv_gas)
    x = ((zz * chit_gas) / (temp * cv_gas));
    //       gam3_gas  = x + 1.0d0
    gam3_gas = (x + 1.0);
    //       gam1_gas  = chit_gas*x + chid_gas
    gam1_gas = ((chit_gas * x) + chid_gas);
    //       nabad_gas = x/gam1_gas
    nabad_gas = (x / gam1_gas);
    //       gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
    gam2_gas = (1.0 / (1.0 - nabad_gas));
    //       cp_gas    = cv_gas * gam1_gas/chid_gas
    cp_gas = ((cv_gas * gam1_gas) / chid_gas);
    //       z         = 1.0d0 + (egas + light2)*zzi
    z = (1.0 + ((egas + light2) * zzi));
    //       sound_gas = clight * sqrt(gam1_gas/z)
    sound_gas = (clight * std::sqrt((gam1_gas / z)));
    //
    //
    //
    //! for the totals
    //       zz    = pres*deni
    zz = (pres * deni);
    //       zzi   = den/pres
    zzi = (den / pres);
    //       chit  = temp/pres * dpresdt
    chit = ((temp / pres) * dpresdt);
    //       chid  = dpresdd*zzi
    chid = (dpresdd * zzi);
    //       cv    = denerdt
    cv = denerdt;
    //       x     = zz * chit/(temp * cv)
    x = ((zz * chit) / (temp * cv));
    //       gam3  = x + 1.0d0
    gam3 = (x + 1.0);
    //       gam1  = chit*x + chid
    gam1 = ((chit * x) + chid);
    //       nabad = x/gam1
    nabad = (x / gam1);
    //       gam2  = 1.0d0/(1.0d0 - nabad)
    gam2 = (1.0 / (1.0 - nabad));
    //       cp    = cv * gam1/chid
    cp = ((cv * gam1) / chid);
    //       z     = 1.0d0 + (ener + light2)*zzi
    z = (1.0 + ((ener + light2) * zzi));
    //       sound = clight * sqrt(gam1/z)
    sound = (clight * std::sqrt((gam1 / z)));
    //
    //
    //
    //! maxwell relations; each is zero if the consistency is perfect
    //       x   = den * den
    x = (den * den);
    //
    //       dse = temp*dentrdt/denerdt - 1.0d0
    dse = (((temp * dentrdt) / denerdt) - 1.0);
    //
    //       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
    dpe = ((((denerdd * x) + (temp * dpresdt)) / pres) - 1.0);
    //
    //       dsp = -dentrdd*x/dpresdt - 1.0d0
    dsp = ((((-dentrdd) * x) / dpresdt) - 1.0);
    //
    //
    //! store this row
    //        ptot_row(j)   = pres
    out.ptot = pres;
    //        dpt_row(j)    = dpresdt
    out.dpt = dpresdt;
    //        dpd_row(j)    = dpresdd
    out.dpd = dpresdd;
    //        dpa_row(j)    = dpresda
    out.dpa = dpresda;
    //        dpz_row(j)    = dpresdz
    out.dpz = dpresdz;
    //
    //        etot_row(j)   = ener
    out.etot = ener;
    //        det_row(j)    = denerdt
    out.det = denerdt;
    //        ded_row(j)    = denerdd
    out.ded = denerdd;
    //        dea_row(j)    = denerda
    out.dea = denerda;
    //        dez_row(j)    = denerdz
    out.dez = denerdz;
    //
    //        stot_row(j)   = entr
    out.stot = entr;
    //        dst_row(j)    = dentrdt
    out.dst = dentrdt;
    //        dsd_row(j)    = dentrdd
    out.dsd = dentrdd;
    //        dsa_row(j)    = dentrda
    out.dsa = dentrda;
    //        dsz_row(j)    = dentrdz
    out.dsz = dentrdz;
    //
    //
    //        pgas_row(j)   = pgas
    out.pgas = pgas;
    //        dpgast_row(j) = dpgasdt
    out.dpgast = dpgasdt;
    //        dpgasd_row(j) = dpgasdd
    out.dpgasd = dpgasdd;
    //        dpgasa_row(j) = dpgasda
    out.dpgasa = dpgasda;
    //        dpgasz_row(j) = dpgasdz
    out.dpgasz = dpgasdz;
    //
    //        egas_row(j)   = egas
    out.egas = egas;
    //        degast_row(j) = degasdt
    out.degast = degasdt;
    //        degasd_row(j) = degasdd
    out.degasd = degasdd;
    //        degasa_row(j) = degasda
    out.degasa = degasda;
    //        degasz_row(j) = degasdz
    out.degasz = degasdz;
    //
    //        sgas_row(j)   = sgas
    out.sgas = sgas;
    //        dsgast_row(j) = dsgasdt
    out.dsgast = dsgasdt;
    //        dsgasd_row(j) = dsgasdd
    out.dsgasd = dsgasdd;
    //        dsgasa_row(j) = dsgasda
    out.dsgasa = dsgasda;
    //        dsgasz_row(j) = dsgasdz
    out.dsgasz = dsgasdz;
    //
    //
    //        prad_row(j)   = prad
    out.prad = prad;
    //        dpradt_row(j) = dpraddt
    out.dpradt = dpraddt;
    //        dpradd_row(j) = dpraddd
    out.dpradd = dpraddd;
    //        dprada_row(j) = dpradda
    out.dprada = dpradda;
    //        dpradz_row(j) = dpraddz
    out.dpradz = dpraddz;
    //
    //        erad_row(j)   = erad
    out.erad = erad;
    //        deradt_row(j) = deraddt
    out.deradt = deraddt;
    //        deradd_row(j) = deraddd
    out.deradd = deraddd;
    //        derada_row(j) = deradda
    out.derada = deradda;
    //        deradz_row(j) = deraddz
    out.deradz = deraddz;
    //
    //        srad_row(j)   = srad
    out.srad = srad;
    //        dsradt_row(j) = dsraddt
    out.dsradt = dsraddt;
    //        dsradd_row(j) = dsraddd
    out.dsradd = dsraddd;
    //        dsrada_row(j) = dsradda
    out.dsrada = dsradda;
    //        dsradz_row(j) = dsraddz
    out.dsradz = dsraddz;
    //
    //
    //        pion_row(j)   = pion
    out.pion = pion;
    //        dpiont_row(j) = dpiondt
    out.dpiont = dpiondt;
    //        dpiond_row(j) = dpiondd
    out.dpiond = dpiondd;
    //        dpiona_row(j) = dpionda
    out.dpiona = dpionda;
    //        dpionz_row(j) = dpiondz
    out.dpionz = dpiondz;
    //
    //        eion_row(j)   = eion
    out.eion = eion;
    //        deiont_row(j) = deiondt
    out.deiont = deiondt;
    //        deiond_row(j) = deiondd
    out.deiond = deiondd;
    //        deiona_row(j) = deionda
    out.deiona = deionda;
    //        deionz_row(j) = deiondz
    out.deionz = deiondz;
    //
    //        sion_row(j)   = sion
    out.sion = sion;
    //        dsiont_row(j) = dsiondt
    out.dsiont = dsiondt;
    //        dsiond_row(j) = dsiondd
    out.dsiond = dsiondd;
    //        dsiona_row(j) = dsionda
    out.dsiona = dsionda;
    //        dsionz_row(j) = dsiondz
    out.dsionz = dsiondz;
    //
    //        xni_row(j)    = xni
    out.xni = xni;
    //
    //        pele_row(j)   = pele
    out.pele = pele;
    //        ppos_row(j)   = 0.0d0
    out.ppos = 0.0;
    //        dpept_row(j)  = dpepdt
    out.dpept = dpepdt;
    //        dpepd_row(j)  = dpepdd
    out.dpepd = dpepdd;
    //        dpepa_row(j)  = dpepda
    out.dpepa = dpepda;
    //        dpepz_row(j)  = dpepdz
    out.dpepz = dpepdz;
    //
    //        eele_row(j)   = eele
    out.eele = eele;
    //        epos_row(j)   = 0.0d0
    out.epos = 0.0;
    //        deept_row(j)  = deepdt
    out.deept = deepdt;
    //        deepd_row(j)  = deepdd
    out.deepd = deepdd;
    //        deepa_row(j)  = deepda
    out.deepa = deepda;
    //        deepz_row(j)  = deepdz
    out.deepz = deepdz;
    //
    //        sele_row(j)   = sele
    out.sele = sele;
    //        spos_row(j)   = 0.0d0
    out.spos = 0.0;
    //        dsept_row(j)  = dsepdt
    out.dsept = dsepdt;
    //        dsepd_row(j)  = dsepdd
    out.dsepd = dsepdd;
    //        dsepa_row(j)  = dsepda
    out.dsepa = dsepda;
    //        dsepz_row(j)  = dsepdz
    out.dsepz = dsepdz;
    //
    //        xnem_row(j)   = xnem
    out.xnem = xnem;
    //        xne_row(j)    = xnefer
    out.xne = xnefer;
    //        dxnet_row(j)  = dxnedt
    out.dxnet = dxnedt;
    //        dxned_row(j)  = dxnedd
    out.dxned = dxnedd;
    //        dxnea_row(j)  = dxneda
    out.dxnea = dxneda;
    //        dxnez_row(j)  = dxnedz
    out.dxnez = dxnedz;
    //        xnp_row(j)    = 0.0d0
    out.xnp = 0.0;
    //        zeff_row(j)   = zbar
    out.zeff = zbar;
    //
    //        etaele_row(j) = etaele
    out.etaele = etaele;
    //        detat_row(j)  = detadt
    out.detat = detadt;
    //        detad_row(j)  = detadd
    out.detad = detadd;
    //        detaa_row(j)  = detada
    out.detaa = detada;
    //        detaz_row(j)  = detadz
    out.detaz = detadz;
    //        etapos_row(j) = 0.0d0
    out.etapos = 0.0;
    //
    //        pcou_row(j)   = pcoul
    out.pcou = pcoul;
    //        dpcout_row(j) = dpcouldt
    out.dpcout = dpcouldt;
    //        dpcoud_row(j) = dpcouldd
    out.dpcoud = dpcouldd;
    //        dpcoua_row(j) = dpcoulda
    out.dpcoua = dpcoulda;
    //        dpcouz_row(j) = dpcouldz
    out.dpcouz = dpcouldz;
    //
    //        ecou_row(j)   = ecoul
    out.ecou = ecoul;
    //        decout_row(j) = decouldt
    out.decout = decouldt;
    //        decoud_row(j) = decouldd
    out.decoud = decouldd;
    //        decoua_row(j) = decoulda
    out.decoua = decoulda;
    //        decouz_row(j) = decouldz
    out.decouz = decouldz;
    //
    //        scou_row(j)   = scoul
    out.scou = scoul;
    //        dscout_row(j) = dscouldt
    out.dscout = dscouldt;
    //        dscoud_row(j) = dscouldd
    out.dscoud = dscouldd;
    //        dscoua_row(j) = dscoulda
    out.dscoua = dscoulda;
    //        dscouz_row(j) = dscouldz
    out.dscouz = dscouldz;
    //
    //        plasg_row(j)  = plasg
    out.plasg = plasg;
    //
    //        dse_row(j)    = dse
    out.dse = dse;
    //        dpe_row(j)    = dpe
    out.dpe = dpe;
    //        dsp_row(j)    = dsp
    out.dsp = dsp;
    //
    //        cv_gas_row(j)    = cv_gas
    out.cv_gas = cv_gas;
    //        cp_gas_row(j)    = cp_gas
    out.cp_gas = cp_gas;
    //        gam1_gas_row(j)  = gam1_gas
    out.gam1_gas = gam1_gas;
    //        gam2_gas_row(j)  = gam2_gas
    out.gam2_gas = gam2_gas;
    //        gam3_gas_row(j)  = gam3_gas
    out.gam3_gas = gam3_gas;
    //        nabad_gas_row(j) = nabad_gas
    out.nabad_gas = nabad_gas;
    //        cs_gas_row(j)    = sound_gas
    out.cs_gas = sound_gas;
    //
    //        cv_row(j)     = cv
    out.cv = cv;
    //        cp_row(j)     = cp
    out.cp = cp;
    //        gam1_row(j)   = gam1
    out.gam1 = gam1;
    //        gam2_row(j)   = gam2
    out.gam2 = gam2;
    //        gam3_row(j)   = gam3
    out.gam3 = gam3;
    //        nabad_row(j)  = nabad
    out.nabad = nabad;
    //        cs_row(j)     = sound
    out.cs = sound;
    //
    //! end of pipeline loop
    //      enddo
    //      return
    //      end
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    return out;
}
} // namespace octotigerII::helmholtz
