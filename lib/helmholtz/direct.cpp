// C++ translation of F. X. (Frank) Timmes's stellar EOS routines.
// Original scientific algorithms: Timmes & Arnett (1999), ApJS 125, 277;
// Timmes & Swesty (2000), ApJS 126, 501. See NOTICE.md for provenance.
// Original Fortran statements are retained immediately above their translations.
#include "direct.hpp"
#include "constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace octotigerII::helmholtz::detail {
namespace {
using Integrand = void (*)(double, const double *, int, double &, double &, double &, double &, double &,
                           double &);
void etages(double xni, double zbar, double temp, double &eta);
void dfermi(double dk, double eta, double theta, double &fd, double &fdeta, double &fdtheta, double &fdeta2,
            double &fdtheta2, double &fdetadtheta);
void fdfunc1(double x, const double *par, [[maybe_unused]] int n, double &fd, double &fdeta, double &fdtheta, double &fdeta2,
             double &fdtheta2, double &fdetadtheta);
void fdfunc2(double x, const double *par, [[maybe_unused]] int n, double &fd, double &fdeta, double &fdtheta, double &fdeta2,
             double &fdtheta2, double &fdetadtheta);
void dqleg020(Integrand f, double a, double b, double &result, double &drdeta, double &drdtheta,
              double &drdeta2, double &drdtheta2, double &drdetadtheta, const double *par, int n);
void dqlag020(Integrand f, double a, double b, double &result, double &drdeta, double &drdtheta,
              double &drdeta2, double &drdtheta2, double &drdetadtheta, const double *par, int n);
struct DirectWork {
    double etaele{};
    double detadd{};
    double detadt{};
    double detada{};
    double detadz{};
    double detaddd{};
    double detaddt{};
    double detadda{};
    double detaddz{};
    double detadtt{};
    double detadta{};
    double detadtz{};
    double detadaa{};
    double detadaz{};
    double detadzz{};
    double pep{};
    double dpepdd{};
    double dpepdt{};
    double dpepda{};
    double dpepdz{};
    double dpepddd{};
    double dpepddt{};
    double dpepdda{};
    double dpepddz{};
    double dpepdtt{};
    double dpepdta{};
    double dpepdtz{};
    double dpepdaa{};
    double dpepdaz{};
    double dpepdzz{};
    double eep{};
    double deepdd{};
    double deepdt{};
    double deepda{};
    double deepdz{};
    double deepddd{};
    double deepddt{};
    double deepdda{};
    double deepddz{};
    double deepdtt{};
    double deepdta{};
    double deepdtz{};
    double deepdaa{};
    double deepdaz{};
    double deepdzz{};
    double sep{};
    double dsepdd{};
    double dsepdt{};
    double dsepda{};
    double dsepdz{};
    double dsepddd{};
    double dsepddt{};
    double dsepdda{};
    double dsepddz{};
    double dsepdtt{};
    double dsepdta{};
    double dsepdtz{};
    double dsepdaa{};
    double dsepdaz{};
    double dsepdzz{};
    double etapos{};
    double zeff{};
    double pele{};
    double dpeledd{};
    double dpeledt{};
    double dpeleda{};
    double dpeledz{};
    double dpeleddd{};
    double dpeleddt{};
    double dpeledda{};
    double dpeleddz{};
    double dpeledtt{};
    double dpeledta{};
    double dpeledtz{};
    double dpeledaa{};
    double dpeledaz{};
    double dpeledzz{};
    double eele{};
    double deeledd{};
    double deeledt{};
    double deeleda{};
    double deeledz{};
    double deeleddd{};
    double deeleddt{};
    double deeledda{};
    double deeleddz{};
    double deeledtt{};
    double deeledta{};
    double deeledtz{};
    double deeledaa{};
    double deeledaz{};
    double deeledzz{};
    double sele{};
    double dseledd{};
    double dseledt{};
    double dseleda{};
    double dseledz{};
    double dseleddd{};
    double dseleddt{};
    double dseledda{};
    double dseleddz{};
    double dseledtt{};
    double dseledta{};
    double dseledtz{};
    double dseledaa{};
    double dseledaz{};
    double dseledzz{};
    double ppos{};
    double dpposdd{};
    double dpposdt{};
    double dpposda{};
    double dpposdz{};
    double dpposddd{};
    double dpposddt{};
    double dpposdda{};
    double dpposddz{};
    double dpposdtt{};
    double dpposdta{};
    double dpposdtz{};
    double dpposdaa{};
    double dpposdaz{};
    double dpposdzz{};
    double epos{};
    double deposdd{};
    double deposdt{};
    double deposda{};
    double deposdz{};
    double deposddd{};
    double deposddt{};
    double deposdda{};
    double deposddz{};
    double deposdtt{};
    double deposdta{};
    double deposdtz{};
    double deposdaa{};
    double deposdaz{};
    double deposdzz{};
    double spos{};
    double dsposdd{};
    double dsposdt{};
    double dsposda{};
    double dsposdz{};
    double dsposddd{};
    double dsposddt{};
    double dsposdda{};
    double dsposddz{};
    double dsposdtt{};
    double dsposdta{};
    double dsposdtz{};
    double dsposdaa{};
    double dsposdaz{};
    double dsposdzz{};
    double xne{};
    double dxnedd{};
    double dxnedt{};
    double dxneda{};
    double dxnedz{};
    double dxneddd{};
    double dxneddt{};
    double dxnedda{};
    double dxneddz{};
    double dxnedtt{};
    double dxnedta{};
    double dxnedtz{};
    double dxnedaa{};
    double dxnedaz{};
    double dxnedzz{};
    double xnefer{};
    double dxneferdd{};
    double dxneferdt{};
    double dxneferda{};
    double dxneferdz{};
    double dxneferddd{};
    double dxneferddt{};
    double dxneferdda{};
    double dxneferddz{};
    double dxneferdtt{};
    double dxneferdta{};
    double dxneferdtz{};
    double dxneferdaa{};
    double dxneferdaz{};
    double dxneferdzz{};
    double xnpfer{};
    double dxnpferdd{};
    double dxnpferdt{};
    double dxnpferda{};
    double dxnpferdz{};
    double dxnpferddd{};
    double dxnpferddt{};
    double dxnpferdda{};
    double dxnpferddz{};
    double dxnpferdtt{};
    double dxnpferdta{};
    double dxnpferdtz{};
    double dxnpferdaa{};
    double dxnpferdaz{};
    double dxnpferdzz{};
    double eip{};
    double deipdd{};
    double deipdt{};
    double deipda{};
    double deipdz{};
    double sip{};
    double dsipdd{};
    double dsipdt{};
    double dsipda{};
    double dsipdz{};
    double pip{};
    double dpipdd{};
    double dpipdt{};
    double dpipda{};
    double dpipdz{};
    void xneroot(int mode, double den, double temp, double abar, double zbar, int ionized, int potmult,
                 double aa, double &f, double &df);
};
void DirectWork::xneroot(int mode, double den, double temp, double abar, double zbar, int ionized,
                         int potmult, double aa, double &f, double &df) {
    [[maybe_unused]] double deni{};
    [[maybe_unused]] double kt{};
    [[maybe_unused]] double kti{};
    [[maybe_unused]] double beta{};
    [[maybe_unused]] double beta12{};
    [[maybe_unused]] double beta32{};
    [[maybe_unused]] double beta52{};
    [[maybe_unused]] double f12{};
    [[maybe_unused]] double f12eta{};
    [[maybe_unused]] double f12beta{};
    [[maybe_unused]] double f12eta2{};
    [[maybe_unused]] double f12beta2{};
    [[maybe_unused]] double f12etabeta{};
    [[maybe_unused]] double f32{};
    [[maybe_unused]] double f32eta{};
    [[maybe_unused]] double f32beta{};
    [[maybe_unused]] double f32eta2{};
    [[maybe_unused]] double f32beta2{};
    [[maybe_unused]] double f32etabeta{};
    [[maybe_unused]] double f52{};
    [[maybe_unused]] double f52eta{};
    [[maybe_unused]] double f52beta{};
    [[maybe_unused]] double f52eta2{};
    [[maybe_unused]] double f52beta2{};
    [[maybe_unused]] double f52etabeta{};
    [[maybe_unused]] double ytot1{};
    [[maybe_unused]] double zz{};
    [[maybe_unused]] double y{};
    [[maybe_unused]] double yy{};
    [[maybe_unused]] double ww{};
    [[maybe_unused]] double dum1{};
    [[maybe_unused]] double dum2{};
    [[maybe_unused]] double dum3{};
    [[maybe_unused]] double denion{};
    [[maybe_unused]] double xni{};
    [[maybe_unused]] double dxnidd{};
    [[maybe_unused]] double dxnidt{};
    [[maybe_unused]] double dxnida{};
    [[maybe_unused]] double dxnidz{};
    [[maybe_unused]] double dxniddd{};
    [[maybe_unused]] double dxniddt{};
    [[maybe_unused]] double dxnidda{};
    [[maybe_unused]] double dxniddz{};
    [[maybe_unused]] double dxnidtt{};
    [[maybe_unused]] double dxnidta{};
    [[maybe_unused]] double dxnidtz{};
    [[maybe_unused]] double dxnidaa{};
    [[maybe_unused]] double dxnidaz{};
    [[maybe_unused]] double dxnidzz{};
    [[maybe_unused]] double chi{};
    [[maybe_unused]] double chifac{};
    [[maybe_unused]] double dchifacdt{};
    [[maybe_unused]] double dchifacdz{};
    [[maybe_unused]] double dchifacdtt{};
    [[maybe_unused]] double dchifacdtz{};
    [[maybe_unused]] double dchifacdzz{};
    [[maybe_unused]] double saha{};
    [[maybe_unused]] double dsaha_dd{};
    [[maybe_unused]] double dsaha_dt{};
    [[maybe_unused]] double dsaha_da{};
    [[maybe_unused]] double dsaha_dz{};
    [[maybe_unused]] double dsaha_deta{};
    [[maybe_unused]] double dsaha_ddd{};
    [[maybe_unused]] double dsaha_ddt{};
    [[maybe_unused]] double dsaha_dda{};
    [[maybe_unused]] double dsaha_ddz{};
    [[maybe_unused]] double dsaha_dtt{};
    [[maybe_unused]] double dsaha_dta{};
    [[maybe_unused]] double dsaha_dtz{};
    [[maybe_unused]] double dsaha_daa{};
    [[maybe_unused]] double dsaha_daz{};
    [[maybe_unused]] double dsaha_dzz{};
    [[maybe_unused]] double dsaha_deta_dd{};
    [[maybe_unused]] double dsaha_deta_dt{};
    [[maybe_unused]] double dsaha_deta_da{};
    [[maybe_unused]] double dsaha_deta_dz{};
    [[maybe_unused]] double dsaha_deta2{};
    [[maybe_unused]] double sfac{};
    [[maybe_unused]] double dsfac_dd{};
    [[maybe_unused]] double dsfac_dt{};
    [[maybe_unused]] double dsfac_da{};
    [[maybe_unused]] double dsfac_dz{};
    [[maybe_unused]] double dsfac_deta{};
    [[maybe_unused]] double dsfac_ddd{};
    [[maybe_unused]] double dsfac_ddt{};
    [[maybe_unused]] double dsfac_dda{};
    [[maybe_unused]] double dsfac_ddz{};
    [[maybe_unused]] double dsfac_dtt{};
    [[maybe_unused]] double dsfac_dta{};
    [[maybe_unused]] double dsfac_dtz{};
    [[maybe_unused]] double dsfac_daa{};
    [[maybe_unused]] double dsfac_daz{};
    [[maybe_unused]] double dsfac_dzz{};
    [[maybe_unused]] double dsfac_deta_dd{};
    [[maybe_unused]] double dsfac_deta_dt{};
    [[maybe_unused]] double dsfac_deta_da{};
    [[maybe_unused]] double dsfac_deta_dz{};
    [[maybe_unused]] double dsfac_deta2{};
    [[maybe_unused]] double dzeff_dd{};
    [[maybe_unused]] double dzeff_dt{};
    [[maybe_unused]] double dzeff_da{};
    [[maybe_unused]] double dzeff_dz{};
    [[maybe_unused]] double dzeff_deta{};
    [[maybe_unused]] double dzeff_ddd{};
    [[maybe_unused]] double dzeff_ddt{};
    [[maybe_unused]] double dzeff_dda{};
    [[maybe_unused]] double dzeff_ddz{};
    [[maybe_unused]] double dzeff_dtt{};
    [[maybe_unused]] double dzeff_dta{};
    [[maybe_unused]] double dzeff_dtz{};
    [[maybe_unused]] double dzeff_daa{};
    [[maybe_unused]] double dzeff_daz{};
    [[maybe_unused]] double dzeff_dzz{};
    [[maybe_unused]] double dzeff_deta_dd{};
    [[maybe_unused]] double dzeff_deta_dt{};
    [[maybe_unused]] double dzeff_deta_da{};
    [[maybe_unused]] double dzeff_deta_dz{};
    [[maybe_unused]] double dzeff_deta2{};
    [[maybe_unused]] double dxne_dd{};
    [[maybe_unused]] double dxne_dt{};
    [[maybe_unused]] double dxne_da{};
    [[maybe_unused]] double dxne_dz{};
    [[maybe_unused]] double dxne_deta{};
    [[maybe_unused]] double dxne_ddd{};
    [[maybe_unused]] double dxne_ddt{};
    [[maybe_unused]] double dxne_dda{};
    [[maybe_unused]] double dxne_ddz{};
    [[maybe_unused]] double dxne_dtt{};
    [[maybe_unused]] double dxne_dta{};
    [[maybe_unused]] double dxne_dtz{};
    [[maybe_unused]] double dxne_daa{};
    [[maybe_unused]] double dxne_daz{};
    [[maybe_unused]] double dxne_dzz{};
    [[maybe_unused]] double dxne_deta_dd{};
    [[maybe_unused]] double dxne_deta_dt{};
    [[maybe_unused]] double dxne_deta_da{};
    [[maybe_unused]] double dxne_deta_dz{};
    [[maybe_unused]] double dxne_deta2{};
    [[maybe_unused]] double dxnefer_deta{};
    [[maybe_unused]] double dxnefer_dbeta{};
    [[maybe_unused]] double dxnefer_deta2{};
    [[maybe_unused]] double dxnefer_dbeta2{};
    [[maybe_unused]] double dxnefer_deta_dbeta{};
    [[maybe_unused]] double dxnpfer_detap{};
    [[maybe_unused]] double dxnpfer_detap2{};
    [[maybe_unused]] double dxnpfer_detap_dbeta{};
    [[maybe_unused]] double detap_deta{};
    [[maybe_unused]] double detap_dbeta{};
    [[maybe_unused]] double detap_deta2{};
    [[maybe_unused]] double detap_dbeta2{};
    [[maybe_unused]] double detap_deta_dbeta{};
    [[maybe_unused]] double dxnpfer_deta{};
    [[maybe_unused]] double dxnpfer_dbeta{};
    [[maybe_unused]] double dxnpfer_deta2{};
    [[maybe_unused]] double dxnpfer_dbeta2{};
    [[maybe_unused]] double dxnpfer_deta_dbeta{};
    [[maybe_unused]] double dxep_deta{};
    [[maybe_unused]] double dxep_dbeta{};
    [[maybe_unused]] double dxep_deta2{};
    [[maybe_unused]] double dxep_dbeta2{};
    [[maybe_unused]] double dxep_deta_dbeta{};
    [[maybe_unused]] double dzeffdd{};
    [[maybe_unused]] double dzeffdt{};
    [[maybe_unused]] double dzeffda{};
    [[maybe_unused]] double dzeffdz{};
    [[maybe_unused]] double dzeffddd{};
    [[maybe_unused]] double dzeffddt{};
    [[maybe_unused]] double dzeffdda{};
    [[maybe_unused]] double dzeffddz{};
    [[maybe_unused]] double dzeffdtt{};
    [[maybe_unused]] double dzeffdta{};
    [[maybe_unused]] double dzeffdtz{};
    [[maybe_unused]] double dzeffdaa{};
    [[maybe_unused]] double dzeffdaz{};
    [[maybe_unused]] double dzeffdzz{};
    [[maybe_unused]] double dpele_deta{};
    [[maybe_unused]] double dpele_dbeta{};
    [[maybe_unused]] double dpele_deta2{};
    [[maybe_unused]] double dpele_dbeta2{};
    [[maybe_unused]] double dpele_deta_dbeta{};
    [[maybe_unused]] double deele_deta{};
    [[maybe_unused]] double deele_dbeta{};
    [[maybe_unused]] double deele_deta2{};
    [[maybe_unused]] double deele_dbeta2{};
    [[maybe_unused]] double deele_deta_dbeta{};
    [[maybe_unused]] double dppos_deta{};
    [[maybe_unused]] double dppos_dbeta{};
    [[maybe_unused]] double dppos_deta2{};
    [[maybe_unused]] double dppos_dbeta2{};
    [[maybe_unused]] double dppos_deta_dbeta{};
    [[maybe_unused]] double dppos_detap{};
    [[maybe_unused]] double dppos_detap2{};
    [[maybe_unused]] double dppos_detap_dbeta{};
    [[maybe_unused]] double depos_deta{};
    [[maybe_unused]] double depos_dbeta{};
    [[maybe_unused]] double depos_deta2{};
    [[maybe_unused]] double depos_dbeta2{};
    [[maybe_unused]] double depos_deta_dbeta{};
    [[maybe_unused]] double depos_detap{};
    [[maybe_unused]] double depos_detap2{};
    [[maybe_unused]] double depos_detap_dbeta{};
    [[maybe_unused]] double deipddd{};
    [[maybe_unused]] double deipddt{};
    [[maybe_unused]] double deipdda{};
    [[maybe_unused]] double deipddz{};
    [[maybe_unused]] double deipdtt{};
    [[maybe_unused]] double deipdta{};
    [[maybe_unused]] double deipdtz{};
    [[maybe_unused]] double deipdaa{};
    [[maybe_unused]] double deipdaz{};
    [[maybe_unused]] double deipdzz{};
    [[maybe_unused]] double dsipddd{};
    [[maybe_unused]] double dsipddt{};
    [[maybe_unused]] double dsipdda{};
    [[maybe_unused]] double dsipddz{};
    [[maybe_unused]] double dsipdtt{};
    [[maybe_unused]] double dsipdta{};
    [[maybe_unused]] double dsipdtz{};
    [[maybe_unused]] double dsipdaa{};
    [[maybe_unused]] double dsipdaz{};
    [[maybe_unused]] double dsipdzz{};
    [[maybe_unused]] double xconst{};
    [[maybe_unused]] double pconst{};
    [[maybe_unused]] double econst{};
    [[maybe_unused]] double mecc{};
    [[maybe_unused]] double dbetadt{};
    [[maybe_unused]] double safe{};
    [[maybe_unused]] double positron_start{};
    //      subroutine xneroot(mode,den,temp,abar,zbar,ionized,potmult,aa, &
    //                         f,df)
    //
    //      include 'implno.dek'
    //      include 'const.dek'
    //
    //! this routine is called by a root finder to find the degeneracy parameter aa
    //! where the number density from a saha equation equals the number density as
    //! computed by the fermi-dirac integrals.
    //
    //! input:
    //! mode    = 0 = root find on electron number density, = 1 = full calculation
    //! temp    = temperature
    //! den     = density
    //! abar    = average weight
    //! zbar    = average charge
    //! aa      = degeneracy parameter (chemical potential/kerg*temp)
    //! ionized = flag to turn off/on ionization contributions
    //! potmult = flag to turn off/on ionization ptential contributions
    //
    //! output
    //! f = value of the function to be zeroed ; xne - xnefer + xnepos = 0
    //! df = derivative with respect to aa of a
    //
    //
    //! declare the pass
    //      integer          mode,ionized,potmult
    //      double precision den,temp,abar,zbar,aa,f,df
    //
    //
    //! bring in common block variables
    //! common block communication with routine xneroot
    //
    //! electron-positrons
    //
    //      double precision etaele,detadd,detadt,detada,detadz, &
    //                       detaddd,detaddt,detadda,detaddz,detadtt, &
    //                       detadta,detadtz,detadaa,detadaz,detadzz
    //
    //      common /xneta/   etaele,detadd,detadt,detada,detadz, &
    //                       detaddd,detaddt,detadda,detaddz,detadtt, &
    //                       detadta,detadtz,detadaa,detadaz,detadzz
    //
    //
    //      double precision &
    //                       pep,dpepdd,dpepdt,dpepda,dpepdz, &
    //                       dpepddd,dpepddt,dpepdda,dpepddz, &
    //                       dpepdtt,dpepdta,dpepdtz,dpepdaa, &
    //                       dpepdaz,dpepdzz, &
    //                       eep,deepdd,deepdt,deepda,deepdz, &
    //                       deepddd,deepddt,deepdda,deepddz, &
    //                       deepdtt,deepdta,deepdtz,deepdaa, &
    //                       deepdaz,deepdzz, &
    //                       sep,dsepdd,dsepdt,dsepda,dsepdz, &
    //                       dsepddd,dsepddt,dsepdda,dsepddz, &
    //                       dsepdtt,dsepdta,dsepdtz,dsepdaa, &
    //                       dsepdaz,dsepdzz
    //
    //      common /epc1/ &
    //                       pep,dpepdd,dpepdt,dpepda,dpepdz, &
    //                       dpepddd,dpepddt,dpepdda,dpepddz, &
    //                       dpepdtt,dpepdta,dpepdtz,dpepdaa, &
    //                       dpepdaz,dpepdzz, &
    //                       eep,deepdd,deepdt,deepda,deepdz, &
    //                       deepddd,deepddt,deepdda,deepddz, &
    //                       deepdtt,deepdta,deepdtz,deepdaa, &
    //                       deepdaz,deepdzz, &
    //                       sep,dsepdd,dsepdt,dsepda,dsepdz, &
    //                       dsepddd,dsepddt,dsepdda,dsepddz, &
    //                       dsepdtt,dsepdta,dsepdtz,dsepdaa, &
    //                       dsepdaz,dsepdzz
    //
    //
    //      double precision &
    //                       etapos,zeff
    //      common /xnec1/ &
    //                       etapos,zeff
    //
    //
    //      double precision &
    //                       pele,dpeledd,dpeledt,dpeleda,dpeledz, &
    //                       dpeleddd,dpeleddt,dpeledda,dpeleddz, &
    //                       dpeledtt,dpeledta,dpeledtz,dpeledaa, &
    //                       dpeledaz,dpeledzz, &
    //                       eele,deeledd,deeledt,deeleda,deeledz, &
    //                       deeleddd,deeleddt,deeledda,deeleddz, &
    //                       deeledtt,deeledta,deeledtz,deeledaa, &
    //                       deeledaz,deeledzz, &
    //                       sele,dseledd,dseledt,dseleda,dseledz, &
    //                       dseleddd,dseleddt,dseledda,dseleddz, &
    //                       dseledtt,dseledta,dseledtz,dseledaa, &
    //                       dseledaz,dseledzz
    //      common /eleth1/ &
    //                       pele,dpeledd,dpeledt,dpeleda,dpeledz, &
    //                       dpeleddd,dpeleddt,dpeledda,dpeleddz, &
    //                       dpeledtt,dpeledta,dpeledtz,dpeledaa, &
    //                       dpeledaz,dpeledzz, &
    //                       eele,deeledd,deeledt,deeleda,deeledz, &
    //                       deeleddd,deeleddt,deeledda,deeleddz, &
    //                       deeledtt,deeledta,deeledtz,deeledaa, &
    //                       deeledaz,deeledzz, &
    //                       sele,dseledd,dseledt,dseleda,dseledz, &
    //                       dseleddd,dseleddt,dseledda,dseleddz, &
    //                       dseledtt,dseledta,dseledtz,dseledaa, &
    //                       dseledaz,dseledzz
    //
    //      double precision &
    //                       ppos,dpposdd,dpposdt,dpposda,dpposdz, &
    //                       dpposddd,dpposddt,dpposdda,dpposddz, &
    //                       dpposdtt,dpposdta,dpposdtz,dpposdaa, &
    //                       dpposdaz,dpposdzz, &
    //                       epos,deposdd,deposdt,deposda,deposdz, &
    //                       deposddd,deposddt,deposdda,deposddz, &
    //                       deposdtt,deposdta,deposdtz,deposdaa, &
    //                       deposdaz,deposdzz, &
    //                       spos,dsposdd,dsposdt,dsposda,dsposdz, &
    //                       dsposddd,dsposddt,dsposdda,dsposddz, &
    //                       dsposdtt,dsposdta,dsposdtz,dsposdaa, &
    //                       dsposdaz,dsposdzz
    //      common /posth1/ &
    //                       ppos,dpposdd,dpposdt,dpposda,dpposdz, &
    //                       dpposddd,dpposddt,dpposdda,dpposddz, &
    //                       dpposdtt,dpposdta,dpposdtz,dpposdaa, &
    //                       dpposdaz,dpposdzz, &
    //                       epos,deposdd,deposdt,deposda,deposdz, &
    //                       deposddd,deposddt,deposdda,deposddz, &
    //                       deposdtt,deposdta,deposdtz,deposdaa, &
    //                       deposdaz,deposdzz, &
    //                       spos,dsposdd,dsposdt,dsposda,dsposdz, &
    //                       dsposddd,dsposddt,dsposdda,dsposddz, &
    //                       dsposdtt,dsposdta,dsposdtz,dsposdaa, &
    //                       dsposdaz,dsposdzz
    //
    //
    //      double precision xne, &
    //                       dxnedd,dxnedt,dxneda,dxnedz, &
    //                       dxneddd,dxneddt,dxnedda,dxneddz, &
    //                       dxnedtt,dxnedta,dxnedtz,dxnedaa, &
    //                       dxnedaz,dxnedzz
    //
    //      common /xnec2/   xne, &
    //                       dxnedd,dxnedt,dxneda,dxnedz, &
    //                       dxneddd,dxneddt,dxnedda,dxneddz, &
    //                       dxnedtt,dxnedta,dxnedtz,dxnedaa, &
    //                       dxnedaz,dxnedzz
    //
    //      double precision &
    //                       xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
    //                       dxneferddd,dxneferddt,dxneferdda,dxneferddz, &
    //                       dxneferdtt,dxneferdta,dxneferdtz,dxneferdaa, &
    //                       dxneferdaz,dxneferdzz, &
    //                       xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
    //                       dxnpferddd,dxnpferddt,dxnpferdda,dxnpferddz, &
    //                       dxnpferdtt,dxnpferdta,dxnpferdtz,dxnpferdaa, &
    //                       dxnpferdaz,dxnpferdzz
    //      common /xnec3/ &
    //                       xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
    //                       dxneferddd,dxneferddt,dxneferdda,dxneferddz, &
    //                       dxneferdtt,dxneferdta,dxneferdtz,dxneferdaa, &
    //                       dxneferdaz,dxneferdzz, &
    //                       xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
    //                       dxnpferddd,dxnpferddt,dxnpferdda,dxnpferddz, &
    //                       dxnpferdtt,dxnpferdta,dxnpferdtz,dxnpferdaa, &
    //                       dxnpferdaz,dxnpferdzz
    //
    //
    //
    //
    //
    //! ionization contributions
    //      double precision eip,deipdd,deipdt,deipda,deipdz, &
    //                       sip,dsipdd,dsipdt,dsipda,dsipdz, &
    //                       pip,dpipdd,dpipdt,dpipda,dpipdz
    //
    //      common /xnec4/   eip,deipdd,deipdt,deipda,deipdz, &
    //                       sip,dsipdd,dsipdt,dsipda,dsipdz, &
    //                       pip,dpipdd,dpipdt,dpipda,dpipdz
    //
    //
    //
    //
    //! local variables
    //      double precision deni,kt,kti,beta,beta12,beta32,beta52, &
    //                       f12,f12eta,f12beta,f12eta2,f12beta2,f12etabeta, &
    //                       f32,f32eta,f32beta,f32eta2,f32beta2,f32etabeta, &
    //                       f52,f52eta,f52beta,f52eta2,f52beta2,f52etabeta, &
    //                       ytot1,zz,y,yy,ww,dum1,dum2,dum3,denion
    //
    //      double precision xni,dxnidd,dxnidt,dxnida,dxnidz, &
    //                       dxniddd,dxniddt,dxnidda,dxniddz, &
    //                       dxnidtt,dxnidta,dxnidtz,dxnidaa, &
    //                       dxnidaz,dxnidzz
    //
    //      double precision chi,chifac,dchifacdt,dchifacdz, &
    //                       dchifacdtt,dchifacdtz,dchifacdzz
    //
    //      double precision saha, &
    //                       dsaha_dd,dsaha_dt,dsaha_da,dsaha_dz,dsaha_deta, &
    //                       dsaha_ddd,dsaha_ddt,dsaha_dda,dsaha_ddz, &
    //                       dsaha_dtt,dsaha_dta,dsaha_dtz,dsaha_daa, &
    //                       dsaha_daz,dsaha_dzz,dsaha_deta_dd,dsaha_deta_dt, &
    //                       dsaha_deta_da,dsaha_deta_dz,dsaha_deta2
    //
    //      double precision sfac, &
    //                       dsfac_dd,dsfac_dt,dsfac_da,dsfac_dz,dsfac_deta, &
    //                       dsfac_ddd,dsfac_ddt,dsfac_dda,dsfac_ddz, &
    //                       dsfac_dtt,dsfac_dta,dsfac_dtz,dsfac_daa, &
    //                       dsfac_daz,dsfac_dzz,dsfac_deta_dd,dsfac_deta_dt, &
    //                       dsfac_deta_da,dsfac_deta_dz,dsfac_deta2
    //
    //      double precision &
    //                       dzeff_dd,dzeff_dt,dzeff_da,dzeff_dz,dzeff_deta, &
    //                       dzeff_ddd,dzeff_ddt,dzeff_dda,dzeff_ddz, &
    //                       dzeff_dtt,dzeff_dta,dzeff_dtz,dzeff_daa, &
    //                       dzeff_daz,dzeff_dzz,dzeff_deta_dd,dzeff_deta_dt, &
    //                       dzeff_deta_da,dzeff_deta_dz,dzeff_deta2
    //
    //      double precision &
    //                       dxne_dd,dxne_dt,dxne_da,dxne_dz,dxne_deta, &
    //                       dxne_ddd,dxne_ddt,dxne_dda,dxne_ddz, &
    //                       dxne_dtt,dxne_dta,dxne_dtz,dxne_daa, &
    //                       dxne_daz,dxne_dzz,dxne_deta_dd,dxne_deta_dt, &
    //                       dxne_deta_da,dxne_deta_dz,dxne_deta2
    //
    //
    //      double precision dxnefer_deta,dxnefer_dbeta, &
    //                       dxnefer_deta2,dxnefer_dbeta2, &
    //                       dxnefer_deta_dbeta, &
    //                       dxnpfer_detap, &
    //                       dxnpfer_detap2, &
    //                       dxnpfer_detap_dbeta, &
    //                       detap_deta,detap_dbeta, &
    //                       detap_deta2,detap_dbeta2, &
    //                       detap_deta_dbeta, &
    //                       dxnpfer_deta,dxnpfer_dbeta, &
    //                       dxnpfer_deta2,dxnpfer_dbeta2, &
    //                       dxnpfer_deta_dbeta
    //
    //      double precision dxep_deta,dxep_dbeta, &
    //                       dxep_deta2,dxep_dbeta2, &
    //                       dxep_deta_dbeta
    //
    //      double precision &
    //                       dzeffdd,dzeffdt,dzeffda,dzeffdz, &
    //                       dzeffddd,dzeffddt,dzeffdda,dzeffddz, &
    //                       dzeffdtt,dzeffdta,dzeffdtz,dzeffdaa, &
    //                       dzeffdaz,dzeffdzz
    //
    //      double precision &
    //                       dpele_deta,dpele_dbeta, &
    //                       dpele_deta2,dpele_dbeta2,dpele_deta_dbeta, &
    //                       deele_deta,deele_dbeta, &
    //                       deele_deta2,deele_dbeta2,deele_deta_dbeta
    //
    //      double precision &
    //                       dppos_deta,dppos_dbeta, &
    //                       dppos_deta2,dppos_dbeta2,dppos_deta_dbeta, &
    //                       dppos_detap,dppos_detap2,dppos_detap_dbeta, &
    //                       depos_deta,depos_dbeta, &
    //                       depos_deta2,depos_dbeta2,depos_deta_dbeta, &
    //                       depos_detap,depos_detap2,depos_detap_dbeta
    //
    //
    //      double precision &
    //                       deipddd,deipddt,deipdda,deipddz, &
    //                       deipdtt,deipdta,deipdtz,deipdaa, &
    //                       deipdaz,deipdzz, &
    //                       dsipddd,dsipddt,dsipdda,dsipddz, &
    //                       dsipdtt,dsipdta,dsipdtz,dsipdaa, &
    //                       dsipdaz,dsipdzz
    //
    //
    //      double precision xconst,pconst,econst,mecc,dbetadt,safe, &
    //                       positron_start
    //      parameter        (xconst  = 8.0d0 * pi * sqrt(2.0d0) * (me/h)**3 * clight**3, &
    //                        pconst  = xconst * 2.0d0/3.0d0 * me * clight**2, &
    //                        econst  = xconst * me * clight**2, &
    //                        mecc    = me * clight * clight, &
    //                        dbetadt = kerg/mecc, &
    //                        safe    = 0.005d0, &
    //                        positron_start = 0.02d0)
    xconst =
        ((((8.0 * pi) * std::sqrt(2.0)) * ((me / h) * (me / h) * (me / h))) * (clight * clight * clight));
    pconst = ((((xconst * 2.0) / 3.0) * me) * (clight * clight));
    econst = ((xconst * me) * (clight * clight));
    mecc = ((me * clight) * clight);
    dbetadt = (kerg / mecc);
    safe = 0.005;
    positron_start = 0.02;
    //
    //! some common factors
    //      ytot1   = 1.0d0/abar
    ytot1 = (1.0 / abar);
    //      deni    = 1.0d0/den
    deni = (1.0 / den);
    //      kt      = kerg * temp
    kt = (kerg * temp);
    //      kti     = 1.0d0/kt
    kti = (1.0 / kt);
    //      beta    = kt/mecc
    beta = (kt / mecc);
    //      beta12  = sqrt(beta)
    beta12 = std::sqrt(beta);
    //      beta32  = beta * beta12
    beta32 = (beta * beta12);
    //      beta52  = beta * beta32
    beta52 = (beta * beta32);
    //      etaele  = aa
    etaele = aa;
    //
    //
    //! ion number density in 1/cm**3
    //      xni     = avo * ytot1 * den
    xni = ((avo * ytot1) * den);
    //      dxnidd  = avo * ytot1
    dxnidd = (avo * ytot1);
    //      dxnidt  = 0.0d0
    dxnidt = 0.0;
    //      dxnida  = -xni * ytot1
    dxnida = ((-xni) * ytot1);
    //      dxnidz  = 0.0d0
    dxnidz = 0.0;
    //
    //      dxniddd = 0.0d0
    dxniddd = 0.0;
    //      dxniddt = 0.0d0
    dxniddt = 0.0;
    //      dxnidda = -dxnidd*ytot1
    dxnidda = ((-dxnidd) * ytot1);
    //      dxniddz = 0.0d0
    dxniddz = 0.0;
    //      dxnidtt = 0.0d0
    dxnidtt = 0.0;
    //      dxnidta = 0.0d0
    dxnidta = 0.0;
    //      dxnidtz = 0.0d0
    dxnidtz = 0.0;
    //      dxnidaa = -2.0d0 * dxnida * ytot1
    dxnidaa = (((-2.0) * dxnida) * ytot1);
    //      dxnidaz = 0.0d0
    dxnidaz = 0.0;
    //      dxnidzz = 0.0d0
    dxnidzz = 0.0;
    //
    //
    //
    //! get the number density of free electrons
    //! saha is the ratio of the ground state to the ionized state
    //! this model is exact for a pure hydrogen composition
    //! denion is a crude pressure ionization model
    //
    //
    //! assume fully ionized
    //      chi        = 0.0d0
    chi = 0.0;
    //      chifac     = 0.0d0
    chifac = 0.0;
    //      dchifacdt  = 0.0d0
    dchifacdt = 0.0;
    //      dchifacdz  = 0.0d0
    dchifacdz = 0.0;
    //      dchifacdtt = 0.0d0
    dchifacdtt = 0.0;
    //      dchifacdtz = 0.0d0
    dchifacdtz = 0.0;
    //      dchifacdzz = 0.0d0
    dchifacdzz = 0.0;
    //      saha       = 0.0d0
    saha = 0.0;
    //      dsaha_dd   = 0.0d0
    dsaha_dd = 0.0;
    //      dsaha_dt   = 0.0d0
    dsaha_dt = 0.0;
    //      dsaha_da   = 0.0d0
    dsaha_da = 0.0;
    //      dsaha_dz   = 0.0d0
    dsaha_dz = 0.0;
    //      dsaha_deta = 0.0d0
    dsaha_deta = 0.0;
    //      dsaha_ddd  = 0.0d0
    dsaha_ddd = 0.0;
    //      dsaha_ddt  = 0.0d0
    dsaha_ddt = 0.0;
    //      dsaha_dda  = 0.0d0
    dsaha_dda = 0.0;
    //      dsaha_ddz  = 0.0d0
    dsaha_ddz = 0.0;
    //      dsaha_dtt  = 0.0d0
    dsaha_dtt = 0.0;
    //      dsaha_dta  = 0.0d0
    dsaha_dta = 0.0;
    //      dsaha_dtz  = 0.0d0
    dsaha_dtz = 0.0;
    //      dsaha_daa  = 0.0d0
    dsaha_daa = 0.0;
    //      dsaha_daz  = 0.0d0
    dsaha_daz = 0.0;
    //      dsaha_dzz  = 0.0d0
    dsaha_dzz = 0.0;
    //      dsaha_deta_dd = 0.0d0
    dsaha_deta_dd = 0.0;
    //      dsaha_deta_dt = 0.0d0
    dsaha_deta_dt = 0.0;
    //      dsaha_deta_da = 0.0d0
    dsaha_deta_da = 0.0;
    //      dsaha_deta_dz = 0.0d0
    dsaha_deta_dz = 0.0;
    //      dsaha_deta2   = 0.0d0
    dsaha_deta2 = 0.0;
    //
    //
    //! do a simple saha approach
    //
    //      if (ionized .eq. 0) then
    if ((ionized == 0)) {
        //
        //       denion     = 0.1d0
        denion = 0.1;
        //       chi        = hion * ev2erg * zbar
        chi = ((hion * ev2erg) * zbar);
        //       chifac     = chi*kti
        chifac = (chi * kti);
        //       dchifacdt  = -chifac/temp
        dchifacdt = ((-chifac) / temp);
        //       dchifacdz  = chifac/zbar
        dchifacdz = (chifac / zbar);
        //       dchifacdtt = -2.0d0*dchifacdt/temp
        dchifacdtt = (((-2.0) * dchifacdt) / temp);
        //       dchifacdtz = -dchifacdz/temp
        dchifacdtz = ((-dchifacdz) / temp);
        //       dchifacdzz = 0.0d0
        dchifacdzz = 0.0;
        //
        //       yy         = chifac - den/denion
        yy = (chifac - (den / denion));
        //
        //! completely neutral; set it for a fake convergence
        //       if (yy .gt. 200.0) then
        if ((yy > 200.0)) {
            //        saha       = 1.0d90
            saha = 1e+90;
            //        f          = 0.0d0
            f = 0.0;
            //        df         = 1.0d0
            df = 1.0;
            //        etaele     = -100.0d0
            etaele = (-100.0);
            //        if (mode .eq. 0) return
            if ((mode == 0))
                return;
            //
            //! ionization possible
            //       else if (yy .gt. -200.0) then
        } else if ((yy > (-200.0))) {
            //        ww   = min(200.0d0,chifac + etaele - den/denion)
            ww = std::min(200.0, ((chifac + etaele) - (den / denion)));
            //        saha = 2.0d0 * exp(ww)
            saha = (2.0 * std::exp(ww));
            //        if (ww .ne. 200.0d0) then
            if ((ww != 200.0)) {
                //
                //! first derivatives
                //         dsaha_dd   = -saha/denion
                dsaha_dd = ((-saha) / denion);
                //         dsaha_dt   = saha * dchifacdt
                dsaha_dt = (saha * dchifacdt);
                //         dsaha_da   = 0.0d0
                dsaha_da = 0.0;
                //         dsaha_dz   = saha * dchifacdz
                dsaha_dz = (saha * dchifacdz);
                //         dsaha_deta = saha
                dsaha_deta = saha;
                //
                //! second derivatives
                //         dsaha_ddd  = -dsaha_dd/denion
                dsaha_ddd = ((-dsaha_dd) / denion);
                //         dsaha_ddt  = -dsaha_dt/denion
                dsaha_ddt = ((-dsaha_dt) / denion);
                //         dsaha_dda  = -dsaha_da/denion
                dsaha_dda = ((-dsaha_da) / denion);
                //         dsaha_ddz  = -dsaha_dz/denion
                dsaha_ddz = ((-dsaha_dz) / denion);
                //         dsaha_dtt  = dsaha_dt*dchifacdt + saha*dchifacdtt
                dsaha_dtt = ((dsaha_dt * dchifacdt) + (saha * dchifacdtt));
                //         dsaha_dta  = dsaha_da*dchifacdt
                dsaha_dta = (dsaha_da * dchifacdt);
                //         dsaha_dtz  = dsaha_dz*dchifacdt + saha*dchifacdtz
                dsaha_dtz = ((dsaha_dz * dchifacdt) + (saha * dchifacdtz));
                //         dsaha_daa  = 0.0d0
                dsaha_daa = 0.0;
                //         dsaha_daz  = 0.0d0
                dsaha_daz = 0.0;
                //         dsaha_dzz  = dsaha_dz*dchifacdz + saha * dchifacdzz
                dsaha_dzz = ((dsaha_dz * dchifacdz) + (saha * dchifacdzz));
                //         dsaha_deta_dd = dsaha_dd
                dsaha_deta_dd = dsaha_dd;
                //         dsaha_deta_dt = dsaha_dt
                dsaha_deta_dt = dsaha_dt;
                //         dsaha_deta_da = dsaha_da
                dsaha_deta_da = dsaha_da;
                //         dsaha_deta_dz = dsaha_dz
                dsaha_deta_dz = dsaha_dz;
                //         dsaha_deta2   = dsaha_deta
                dsaha_deta2 = dsaha_deta;
                //
                //        end if
            }
            //       end if
        }
        //      end if
    }
    //
    //
    //
    //! saha factor, the fraction ionized
    //      sfac       = 1.0d0/(1.0d0 + saha)
    sfac = (1.0 / (1.0 + saha));
    //
    //! first derivatives
    //      y          = -sfac*sfac
    y = ((-sfac) * sfac);
    //      dsfac_dd   = y*dsaha_dd
    dsfac_dd = (y * dsaha_dd);
    //      dsfac_dt   = y*dsaha_dt
    dsfac_dt = (y * dsaha_dt);
    //      dsfac_da   = y*dsaha_da
    dsfac_da = (y * dsaha_da);
    //      dsfac_dz   = y*dsaha_dz
    dsfac_dz = (y * dsaha_dz);
    //      dsfac_deta = y*dsaha_deta
    dsfac_deta = (y * dsaha_deta);
    //
    //! second derivatives
    //      ww         = -2.0d0*sfac
    ww = ((-2.0) * sfac);
    //      dsfac_ddd  = ww*dsfac_dd*dsaha_dd + y*dsaha_ddd
    dsfac_ddd = (((ww * dsfac_dd) * dsaha_dd) + (y * dsaha_ddd));
    //      dsfac_ddt  = ww*dsfac_dt*dsaha_dd + y*dsaha_ddt
    dsfac_ddt = (((ww * dsfac_dt) * dsaha_dd) + (y * dsaha_ddt));
    //      dsfac_dda  = ww*dsfac_da*dsaha_dd + y*dsaha_dda
    dsfac_dda = (((ww * dsfac_da) * dsaha_dd) + (y * dsaha_dda));
    //      dsfac_ddz  = ww*dsfac_dz*dsaha_dd + y*dsaha_ddz
    dsfac_ddz = (((ww * dsfac_dz) * dsaha_dd) + (y * dsaha_ddz));
    //      dsfac_dtt  = ww*dsfac_dt*dsaha_dt + y*dsaha_dtt
    dsfac_dtt = (((ww * dsfac_dt) * dsaha_dt) + (y * dsaha_dtt));
    //      dsfac_dta  = ww*dsfac_da*dsaha_dt + y*dsaha_dta
    dsfac_dta = (((ww * dsfac_da) * dsaha_dt) + (y * dsaha_dta));
    //      dsfac_dtz  = ww*dsfac_dz*dsaha_dt + y*dsaha_dtz
    dsfac_dtz = (((ww * dsfac_dz) * dsaha_dt) + (y * dsaha_dtz));
    //      dsfac_daa  = ww*dsfac_da*dsaha_da + y*dsaha_daa
    dsfac_daa = (((ww * dsfac_da) * dsaha_da) + (y * dsaha_daa));
    //      dsfac_daz  = ww*dsfac_dz*dsaha_da + y*dsaha_daz
    dsfac_daz = (((ww * dsfac_dz) * dsaha_da) + (y * dsaha_daz));
    //      dsfac_dzz  = ww*dsfac_dz*dsaha_dz + y*dsaha_dzz
    dsfac_dzz = (((ww * dsfac_dz) * dsaha_dz) + (y * dsaha_dzz));
    //      dsfac_deta_dd = ww*dsfac_dd*dsaha_deta + y*dsaha_deta_dd
    dsfac_deta_dd = (((ww * dsfac_dd) * dsaha_deta) + (y * dsaha_deta_dd));
    //      dsfac_deta_dt = ww*dsfac_dt*dsaha_deta + y*dsaha_deta_dt
    dsfac_deta_dt = (((ww * dsfac_dt) * dsaha_deta) + (y * dsaha_deta_dt));
    //      dsfac_deta_da = ww*dsfac_da*dsaha_deta + y*dsaha_deta_da
    dsfac_deta_da = (((ww * dsfac_da) * dsaha_deta) + (y * dsaha_deta_da));
    //      dsfac_deta_dz = ww*dsfac_dz*dsaha_deta + y*dsaha_deta_dz
    dsfac_deta_dz = (((ww * dsfac_dz) * dsaha_deta) + (y * dsaha_deta_dz));
    //      dsfac_deta2   = ww*dsfac_deta*dsaha_deta + y*dsaha_deta2
    dsfac_deta2 = (((ww * dsfac_deta) * dsaha_deta) + (y * dsaha_deta2));
    //
    //
    //
    //! effective charge
    //      zeff       = zbar * sfac
    zeff = (zbar * sfac);
    //
    //! first derivatives
    //      dzeff_dd   = zbar * dsfac_dd
    dzeff_dd = (zbar * dsfac_dd);
    //      dzeff_dt   = zbar * dsfac_dt
    dzeff_dt = (zbar * dsfac_dt);
    //      dzeff_da   = zbar * dsfac_da
    dzeff_da = (zbar * dsfac_da);
    //      dzeff_dz   = sfac
    dzeff_dz = sfac;
    //      dzeff_deta = zbar * dsfac_deta
    dzeff_deta = (zbar * dsfac_deta);
    //
    //! second derivatives
    //      dzeff_ddd  = zbar*dsfac_ddd
    dzeff_ddd = (zbar * dsfac_ddd);
    //      dzeff_ddt  = zbar*dsfac_ddt
    dzeff_ddt = (zbar * dsfac_ddt);
    //      dzeff_dda  = zbar*dsfac_dda
    dzeff_dda = (zbar * dsfac_dda);
    //      dzeff_ddz  = zbar*dsfac_ddz
    dzeff_ddz = (zbar * dsfac_ddz);
    //      dzeff_dtt  = zbar*dsfac_dtt
    dzeff_dtt = (zbar * dsfac_dtt);
    //      dzeff_dta  = zbar*dsfac_dta
    dzeff_dta = (zbar * dsfac_dta);
    //      dzeff_dtz  = zbar*dsfac_dtz
    dzeff_dtz = (zbar * dsfac_dtz);
    //      dzeff_daa  = zbar*dsfac_daa
    dzeff_daa = (zbar * dsfac_daa);
    //      dzeff_daz  = zbar*dsfac_daz
    dzeff_daz = (zbar * dsfac_daz);
    //      dzeff_dzz  = dsfac_dz
    dzeff_dzz = dsfac_dz;
    //      dzeff_deta_dd = zbar*dsfac_deta_dd
    dzeff_deta_dd = (zbar * dsfac_deta_dd);
    //      dzeff_deta_dt = zbar*dsfac_deta_dt
    dzeff_deta_dt = (zbar * dsfac_deta_dt);
    //      dzeff_deta_da = zbar*dsfac_deta_da
    dzeff_deta_da = (zbar * dsfac_deta_da);
    //      dzeff_deta_dz = dsfac_deta + zbar*dsfac_deta_dz
    dzeff_deta_dz = (dsfac_deta + (zbar * dsfac_deta_dz));
    //      dzeff_deta2   = zbar * dsfac_deta2
    dzeff_deta2 = (zbar * dsfac_deta2);
    //
    //
    //
    //
    //! number density of free electrons
    //      xne        = xni * zeff
    xne = (xni * zeff);
    //
    //! first derivatives
    //      dxne_dd    = dxnidd * zeff + xni*dzeff_dd
    dxne_dd = ((dxnidd * zeff) + (xni * dzeff_dd));
    //      dxne_dt    = dxnidt * zeff + xni*dzeff_dt
    dxne_dt = ((dxnidt * zeff) + (xni * dzeff_dt));
    //      dxne_da    = dxnida * zeff + xni*dzeff_da
    dxne_da = ((dxnida * zeff) + (xni * dzeff_da));
    //      dxne_dz    = dxnidz * zeff + xni*dzeff_dz
    dxne_dz = ((dxnidz * zeff) + (xni * dzeff_dz));
    //      dxne_deta  = xni*dzeff_deta
    dxne_deta = (xni * dzeff_deta);
    //
    //! second derivatives
    //      dxne_ddd = dxniddd*zeff + 2.0d0*dxnidd*dzeff_dd + xni*dzeff_ddd
    dxne_ddd = (((dxniddd * zeff) + ((2.0 * dxnidd) * dzeff_dd)) + (xni * dzeff_ddd));
    //      dxne_ddt = dxniddt*zeff + dxnidd*dzeff_dt &
    //                 + dxnidt*dzeff_dt + xni*dzeff_ddt
    dxne_ddt = ((((dxniddt * zeff) + (dxnidd * dzeff_dt)) + (dxnidt * dzeff_dt)) + (xni * dzeff_ddt));
    //      dxne_dda = dxnidda*zeff + dxnidd*dzeff_da &
    //                 + dxnida*dzeff_dt + xni*dzeff_dda
    dxne_dda = ((((dxnidda * zeff) + (dxnidd * dzeff_da)) + (dxnida * dzeff_dt)) + (xni * dzeff_dda));
    //      dxne_ddz = dxniddz*zeff + dxnidd*dzeff_dz &
    //                 + dxnidz*dzeff_dt + xni*dzeff_ddz
    dxne_ddz = ((((dxniddz * zeff) + (dxnidd * dzeff_dz)) + (dxnidz * dzeff_dt)) + (xni * dzeff_ddz));
    //      dxne_dtt = dxnidtt*zeff + 2.0d0*dxnidt*dzeff_dt + xni*dzeff_dtt
    dxne_dtt = (((dxnidtt * zeff) + ((2.0 * dxnidt) * dzeff_dt)) + (xni * dzeff_dtt));
    //      dxne_dta = dxnidta*zeff + dxnidt*dzeff_da &
    //                 + dxnida*dzeff_dt + xni*dzeff_dta
    dxne_dta = ((((dxnidta * zeff) + (dxnidt * dzeff_da)) + (dxnida * dzeff_dt)) + (xni * dzeff_dta));
    //      dxne_dtz = dxnidtz*zeff + dxnidt*dzeff_dz &
    //                 + dxnidz*dzeff_dt + xni*dzeff_dtz
    dxne_dtz = ((((dxnidtz * zeff) + (dxnidt * dzeff_dz)) + (dxnidz * dzeff_dt)) + (xni * dzeff_dtz));
    //      dxne_daa = dxnidaa*zeff + 2.0d0*dxnida*dzeff_da + xni*dzeff_daa
    dxne_daa = (((dxnidaa * zeff) + ((2.0 * dxnida) * dzeff_da)) + (xni * dzeff_daa));
    //      dxne_daz = dxnidaz*zeff + dxnida*dzeff_dz &
    //                 + dxnidz*dzeff_da + xni*dzeff_daz
    dxne_daz = ((((dxnidaz * zeff) + (dxnida * dzeff_dz)) + (dxnidz * dzeff_da)) + (xni * dzeff_daz));
    //      dxne_dzz = dxnidzz*zeff + 2.0d0*dxnidz*dzeff_dz + xni*dzeff_dzz
    dxne_dzz = (((dxnidzz * zeff) + ((2.0 * dxnidz) * dzeff_dz)) + (xni * dzeff_dzz));
    //      dxne_deta_dd = dxnidd*dzeff_deta + xni*dzeff_deta_dd
    dxne_deta_dd = ((dxnidd * dzeff_deta) + (xni * dzeff_deta_dd));
    //      dxne_deta_dt = dxnidt*dzeff_deta + xni*dzeff_deta_dt
    dxne_deta_dt = ((dxnidt * dzeff_deta) + (xni * dzeff_deta_dt));
    //      dxne_deta_dz = dxnidz*dzeff_deta + xni*dzeff_deta_dz
    dxne_deta_dz = ((dxnidz * dzeff_deta) + (xni * dzeff_deta_dz));
    //      dxne_deta_da = dxnida*dzeff_deta + xni*dzeff_deta_da
    dxne_deta_da = ((dxnida * dzeff_deta) + (xni * dzeff_deta_da));
    //      dxne_deta2   = xni*dzeff_deta2
    dxne_deta2 = (xni * dzeff_deta2);
    //
    //
    //
    //! get the fermi-dirac integral electron contribution
    //      call dfermi(0.5d0, etaele, beta, f12, f12eta, f12beta, &
    //                                       f12eta2, f12beta2, f12etabeta)
    dfermi(0.5, etaele, beta, f12, f12eta, f12beta, f12eta2, f12beta2, f12etabeta);
    //      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta, &
    //                                       f32eta2, f32beta2, f32etabeta)
    dfermi(1.5, etaele, beta, f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta);
    //
    //      zz   = xconst * beta32
    zz = (xconst * beta32);
    //      ww   = xconst * 1.5d0 * beta12
    ww = ((xconst * 1.5) * beta12);
    //      yy   = f12 + beta * f32
    yy = (f12 + (beta * f32));
    //      dum1 = f12eta + beta*f32eta
    dum1 = (f12eta + (beta * f32eta));
    //      dum2 = f12beta + f32 + beta*f32beta
    dum2 = ((f12beta + f32) + (beta * f32beta));
    //
    //      xnefer             = zz * yy
    xnefer = (zz * yy);
    //      dxnefer_deta       = zz * dum1
    dxnefer_deta = (zz * dum1);
    //      dxnefer_dbeta      = ww*yy + zz*dum2
    dxnefer_dbeta = ((ww * yy) + (zz * dum2));
    //      dxnefer_deta2      = zz * (f12eta2 + beta * f32eta2)
    dxnefer_deta2 = (zz * (f12eta2 + (beta * f32eta2)));
    //      dxnefer_dbeta2     = 0.5d0*ww/beta*yy + 2.0d0*ww*dum2 &
    //                         + zz*(f12beta2 + 2.0d0*f32beta + beta*f32beta2)
    dxnefer_dbeta2 = (((((0.5 * ww) / beta) * yy) + ((2.0 * ww) * dum2)) +
                      (zz * ((f12beta2 + (2.0 * f32beta)) + (beta * f32beta2))));
    //      dxnefer_deta_dbeta = ww*dum1 &
    //                         + zz*(f12etabeta + f32eta + beta*f32etabeta)
    dxnefer_deta_dbeta = ((ww * dum1) + (zz * ((f12etabeta + f32eta) + (beta * f32etabeta))));
    //
    //
    //
    //! if the temperature is not too low, get the positron contributions
    //! chemical equilibrium means etaele + etapos = eta_photon = 0.
    //      etapos             = 0.0d0
    etapos = 0.0;
    //      detap_deta         = 0.0d0
    detap_deta = 0.0;
    //      detap_dbeta        = 0.0d0
    detap_dbeta = 0.0;
    //      detap_deta2        = 0.0d0
    detap_deta2 = 0.0;
    //      detap_dbeta2       = 0.0d0
    detap_dbeta2 = 0.0;
    //      detap_deta_dbeta   = 0.0d0
    detap_deta_dbeta = 0.0;
    //      xnpfer             = 0.0d0
    xnpfer = 0.0;
    //      dxnpfer_deta       = 0.0d0
    dxnpfer_deta = 0.0;
    //      dxnpfer_dbeta      = 0.0d0
    dxnpfer_dbeta = 0.0;
    //      dxnpfer_deta2      = 0.0d0
    dxnpfer_deta2 = 0.0;
    //      dxnpfer_dbeta2     = 0.0d0
    dxnpfer_dbeta2 = 0.0;
    //      dxnpfer_deta_dbeta = 0.0d0
    dxnpfer_deta_dbeta = 0.0;
    //
    //      if (beta .gt. positron_start) then
    if ((beta > positron_start)) {
        //       etapos           = -aa - 2.0d0/beta
        etapos = ((-aa) - (2.0 / beta));
        //       detap_deta       = -1.0d0
        detap_deta = (-1.0);
        //       detap_deta2      = 0.0d0
        detap_deta2 = 0.0;
        //       detap_dbeta      = 2.0d0/beta**2
        detap_dbeta = (2.0 / (beta * beta));
        //       detap_dbeta2     = -4.0d0/beta**3
        detap_dbeta2 = ((-4.0) / (beta * beta * beta));
        //       detap_deta_dbeta = 0.0d0
        detap_deta_dbeta = 0.0;
        //
        //       call dfermi(0.5d0, etapos, beta, f12, f12eta, f12beta, &
        //                                        f12eta2, f12beta2, f12etabeta)
        dfermi(0.5, etapos, beta, f12, f12eta, f12beta, f12eta2, f12beta2, f12etabeta);
        //       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta, &
        //                                        f32eta2, f32beta2, f32etabeta)
        dfermi(1.5, etapos, beta, f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta);
        //
        //       zz   = xconst * beta32
        zz = (xconst * beta32);
        //       ww   = xconst * 1.5d0 * beta12
        ww = ((xconst * 1.5) * beta12);
        //       yy   = f12 + beta * f32
        yy = (f12 + (beta * f32));
        //       dum1 = f12eta + beta*f32eta
        dum1 = (f12eta + (beta * f32eta));
        //       dum2 = f12beta + f32 + beta*f32beta
        dum2 = ((f12beta + f32) + (beta * f32beta));
        //
        //       xnpfer              = zz * yy
        xnpfer = (zz * yy);
        //       dxnpfer_detap       = zz * dum1
        dxnpfer_detap = (zz * dum1);
        //       dxnpfer_dbeta      = ww*yy + zz*dum2
        dxnpfer_dbeta = ((ww * yy) + (zz * dum2));
        //       dxnpfer_detap2      = zz * (f12eta2 + beta * f32eta2)
        dxnpfer_detap2 = (zz * (f12eta2 + (beta * f32eta2)));
        //       dxnpfer_dbeta2     = 0.5d0*ww/beta*yy + 2.0d0*ww*dum2 &
        //                         + zz*(f12beta2 + 2.0d0*f32beta + beta*f32beta2)
        dxnpfer_dbeta2 = (((((0.5 * ww) / beta) * yy) + ((2.0 * ww) * dum2)) +
                          (zz * ((f12beta2 + (2.0 * f32beta)) + (beta * f32beta2))));
        //       dxnpfer_detap_dbeta = ww*dum1 &
        //                         + zz*(f12etabeta + f32eta + beta*f32etabeta)
        dxnpfer_detap_dbeta = ((ww * dum1) + (zz * ((f12etabeta + f32eta) + (beta * f32etabeta))));
        //
        //
        //! convert the etap derivatives to eta derivatives
        //! all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta
        //
        //       dxnpfer_deta  = dxnpfer_detap * detap_deta
        dxnpfer_deta = (dxnpfer_detap * detap_deta);
        //       dxnpfer_dbeta = dxnpfer_dbeta + dxnpfer_detap * detap_dbeta
        dxnpfer_dbeta = (dxnpfer_dbeta + (dxnpfer_detap * detap_dbeta));
        //       dxnpfer_deta2 = dxnpfer_detap2 * detap_deta**2 &
        //                       + dxnpfer_detap * detap_deta2
        dxnpfer_deta2 = ((dxnpfer_detap2 * (detap_deta * detap_deta)) + (dxnpfer_detap * detap_deta2));
        //       dxnpfer_dbeta2 = dxnpfer_dbeta2 &
        //                       + 2.0d0 * dxnpfer_detap_dbeta * detap_dbeta &
        //                       + dxnpfer_detap2 * detap_dbeta**2 &
        //                       + dxnpfer_detap * detap_dbeta2
        dxnpfer_dbeta2 = (((dxnpfer_dbeta2 + ((2.0 * dxnpfer_detap_dbeta) * detap_dbeta)) +
                           (dxnpfer_detap2 * (detap_dbeta * detap_dbeta))) +
                          (dxnpfer_detap * detap_dbeta2));
        //       dxnpfer_deta_dbeta = dxnpfer_detap2 * detap_dbeta * detap_deta &
        //                       + dxnpfer_detap_dbeta * detap_deta &
        //                       + dxnpfer_detap * detap_deta_dbeta
        dxnpfer_deta_dbeta =
            ((((dxnpfer_detap2 * detap_dbeta) * detap_deta) + (dxnpfer_detap_dbeta * detap_deta)) +
             (dxnpfer_detap * detap_deta_dbeta));
        //
        //      end if
    }
    //
    //
    //
    //! charge neutrality means ne_ionization = ne_electrons - ne_positrons
    //      f  = xnefer - xnpfer - xne
    f = ((xnefer - xnpfer) - xne);
    //
    //
    //! derivative of f with eta for newton-like root finders
    //      df = dxnefer_deta  - dxnpfer_deta  - dxne_deta
    df = ((dxnefer_deta - dxnpfer_deta) - dxne_deta);
    //
    //
    //
    //! if we are in root finder mode, return
    //
    //
    //      if (mode .eq. 0) return
    if ((mode == 0))
        return;
    //
    //
    //
    //
    //
    //! if we are not in root finder mode, polish off the calculation
    //
    //
    //! all the derivatives are in terms of eta and beta.
    //! we want to convert to temperature, density, abar and zbar derivatives.
    //! so, after the root find above on eta we have xne = xnefer - xnpfer.
    //! taking the derivative of this for property p
    //! dxne/deta deta + dxne/dp dp = dxnefer/deta deta + dxnefer/dbeta dbeta
    //! solving for the unknown eta derivative yields
    //! deta/dp = (dxne_dp - dxnefer/dbeta dbeta/dp) / (dxnefer/deta - dxne/deta)
    //
    //
    //      dxep_deta       = dxnefer_deta  - dxnpfer_deta
    dxep_deta = (dxnefer_deta - dxnpfer_deta);
    //      dxep_dbeta      = dxnefer_dbeta - dxnpfer_dbeta
    dxep_dbeta = (dxnefer_dbeta - dxnpfer_dbeta);
    //      dxep_deta2      = dxnefer_deta2  - dxnpfer_deta2
    dxep_deta2 = (dxnefer_deta2 - dxnpfer_deta2);
    //      dxep_dbeta2     = dxnefer_dbeta2 - dxnpfer_dbeta2
    dxep_dbeta2 = (dxnefer_dbeta2 - dxnpfer_dbeta2);
    //      dxep_deta_dbeta = dxnefer_deta_dbeta  - dxnpfer_deta_dbeta
    dxep_deta_dbeta = (dxnefer_deta_dbeta - dxnpfer_deta_dbeta);
    //
    //      y          = 1.0d0/(dxep_deta - dxne_deta)
    y = (1.0 / (dxep_deta - dxne_deta));
    //
    //! the all important first derivatives of eta
    //      detadd = dxne_dd * y
    detadd = (dxne_dd * y);
    //      detadt = (dxne_dt - dxep_dbeta*dbetadt) * y
    detadt = ((dxne_dt - (dxep_dbeta * dbetadt)) * y);
    //      detada = dxne_da * y
    detada = (dxne_da * y);
    //      detadz = dxne_dz * y
    detadz = (dxne_dz * y);
    //
    //! second derivatives
    //      detaddd = dxne_ddd*y - detadd*y*(dxep_deta2 - dxne_deta2)*detadd
    detaddd = ((dxne_ddd * y) - (((detadd * y) * (dxep_deta2 - dxne_deta2)) * detadd));
    //      detaddt = dxne_ddt*y - detadd*y*(dxep_deta2*detadt &
    //                + dxep_deta_dbeta*dbetadt - dxne_deta2*detadt)
    detaddt =
        ((dxne_ddt * y) -
         ((detadd * y) * (((dxep_deta2 * detadt) + (dxep_deta_dbeta * dbetadt)) - (dxne_deta2 * detadt))));
    //      detadda = dxne_dda*y - detadd*y*(dxep_deta2 - dxne_deta2)*detada
    detadda = ((dxne_dda * y) - (((detadd * y) * (dxep_deta2 - dxne_deta2)) * detada));
    //      detaddz = dxne_ddz*y - detadd*y*(dxep_deta2 - dxne_deta2)*detadz
    detaddz = ((dxne_ddz * y) - (((detadd * y) * (dxep_deta2 - dxne_deta2)) * detadz));
    //
    //      detadtt = (dxne_dtt - (dxep_deta_dbeta*detadt &
    //                 + dxep_dbeta2*dbetadt)*dbetadt)*y &
    //                 - detadt*y*(dxep_deta2*detadt &
    //                         + dxep_deta_dbeta*dbetadt - dxne_deta2*detadt)
    detadtt =
        (((dxne_dtt - (((dxep_deta_dbeta * detadt) + (dxep_dbeta2 * dbetadt)) * dbetadt)) * y) -
         ((detadt * y) * (((dxep_deta2 * detadt) + (dxep_deta_dbeta * dbetadt)) - (dxne_deta2 * detadt))));
    //
    //      detadta = (dxne_dta - dxep_deta_dbeta*detada*dbetadt)*y &
    //                 - detadt*y*(dxep_deta2*detada - dxne_deta2*detada)
    detadta = (((dxne_dta - ((dxep_deta_dbeta * detada) * dbetadt)) * y) -
               ((detadt * y) * ((dxep_deta2 * detada) - (dxne_deta2 * detada))));
    //
    //      detadtz = (dxne_dtz - dxep_deta_dbeta*detadz*dbetadt)*y &
    //                 - detadt*y*(dxep_deta2*detadz - dxne_deta2*detadz)
    detadtz = (((dxne_dtz - ((dxep_deta_dbeta * detadz) * dbetadt)) * y) -
               ((detadt * y) * ((dxep_deta2 * detadz) - (dxne_deta2 * detadz))));
    //
    //      detadaa = dxne_daa*y - detada*y*(dxep_deta2 - dxne_deta2)*detada
    detadaa = ((dxne_daa * y) - (((detada * y) * (dxep_deta2 - dxne_deta2)) * detada));
    //      detadaz = dxne_daz*y - detada*y*(dxep_deta2 - dxne_deta2)*detadz
    detadaz = ((dxne_daz * y) - (((detada * y) * (dxep_deta2 - dxne_deta2)) * detadz));
    //      detadzz = dxne_dzz*y - detadz*y*(dxep_deta2 - dxne_deta2)*detadz
    detadzz = ((dxne_dzz * y) - (((detadz * y) * (dxep_deta2 - dxne_deta2)) * detadz));
    //
    //
    //
    //! first derivatives of the effective charge
    //      dzeffdd = dzeff_deta*detadd + dzeff_dd
    dzeffdd = ((dzeff_deta * detadd) + dzeff_dd);
    //      dzeffdt = dzeff_deta*detadt + dzeff_dt
    dzeffdt = ((dzeff_deta * detadt) + dzeff_dt);
    //      dzeffda = dzeff_deta*detada + dzeff_da
    dzeffda = ((dzeff_deta * detada) + dzeff_da);
    //      dzeffdz = dzeff_deta*detadz + dzeff_dz
    dzeffdz = ((dzeff_deta * detadz) + dzeff_dz);
    //
    //! second derivatives
    //      dzeffddd = dzeff_deta_dd*detadd + dzeff_deta*detaddd + dzeff_ddd
    dzeffddd = (((dzeff_deta_dd * detadd) + (dzeff_deta * detaddd)) + dzeff_ddd);
    //      dzeffddt = dzeff_deta_dt*detadd + dzeff_deta*detaddt + dzeff_ddt
    dzeffddt = (((dzeff_deta_dt * detadd) + (dzeff_deta * detaddt)) + dzeff_ddt);
    //      dzeffdda = dzeff_deta_da*detadd + dzeff_deta*detadda + dzeff_dda
    dzeffdda = (((dzeff_deta_da * detadd) + (dzeff_deta * detadda)) + dzeff_dda);
    //      dzeffddz = dzeff_deta_dz*detadd + dzeff_deta*detaddz + dzeff_ddz
    dzeffddz = (((dzeff_deta_dz * detadd) + (dzeff_deta * detaddz)) + dzeff_ddz);
    //      dzeffdtt = dzeff_deta_dt*detadt + dzeff_deta*detadtt + dzeff_dtt
    dzeffdtt = (((dzeff_deta_dt * detadt) + (dzeff_deta * detadtt)) + dzeff_dtt);
    //      dzeffdta = dzeff_deta_da*detadt + dzeff_deta*detadta + dzeff_dta
    dzeffdta = (((dzeff_deta_da * detadt) + (dzeff_deta * detadta)) + dzeff_dta);
    //      dzeffdtz = dzeff_deta_dz*detadt + dzeff_deta*detadtz + dzeff_dtz
    dzeffdtz = (((dzeff_deta_dz * detadt) + (dzeff_deta * detadtz)) + dzeff_dtz);
    //      dzeffdaa = dzeff_deta_da*detada + dzeff_deta*detadaa + dzeff_daa
    dzeffdaa = (((dzeff_deta_da * detada) + (dzeff_deta * detadaa)) + dzeff_daa);
    //      dzeffdaz = dzeff_deta_dz*detada + dzeff_deta*detadaz + dzeff_daz
    dzeffdaz = (((dzeff_deta_dz * detada) + (dzeff_deta * detadaz)) + dzeff_daz);
    //      dzeffdzz = dzeff_deta_dz*detadz + dzeff_deta*detadzz + dzeff_dzz
    dzeffdzz = (((dzeff_deta_dz * detadz) + (dzeff_deta * detadzz)) + dzeff_dzz);
    //
    //
    //! first derivatives of the electron number density
    //      dxnedd = dxnidd * zeff + xni * dzeffdd
    dxnedd = ((dxnidd * zeff) + (xni * dzeffdd));
    //      dxnedt = dxnidt * zeff + xni * dzeffdt
    dxnedt = ((dxnidt * zeff) + (xni * dzeffdt));
    //      dxneda = dxnida * zeff + xni * dzeffda
    dxneda = ((dxnida * zeff) + (xni * dzeffda));
    //      dxnedz = dxnidz * zeff + xni * dzeffdz
    dxnedz = ((dxnidz * zeff) + (xni * dzeffdz));
    //
    //! second derivatives
    //      dxneddd = dxniddd*zeff+dxnidd*dzeffdd+dxnidd*dzeffdd+xni*dzeffddd
    dxneddd = ((((dxniddd * zeff) + (dxnidd * dzeffdd)) + (dxnidd * dzeffdd)) + (xni * dzeffddd));
    //      dxneddt = dxniddt*zeff+dxnidd*dzeffdt+dxnidt*dzeffdd+xni*dzeffddt
    dxneddt = ((((dxniddt * zeff) + (dxnidd * dzeffdt)) + (dxnidt * dzeffdd)) + (xni * dzeffddt));
    //      dxnedda = dxnidda*zeff+dxnidd*dzeffda+dxnida*dzeffdd+xni*dzeffdda
    dxnedda = ((((dxnidda * zeff) + (dxnidd * dzeffda)) + (dxnida * dzeffdd)) + (xni * dzeffdda));
    //      dxneddz = dxniddz*zeff+dxnidd*dzeffdz+dxnidz*dzeffdd+xni*dzeffddz
    dxneddz = ((((dxniddz * zeff) + (dxnidd * dzeffdz)) + (dxnidz * dzeffdd)) + (xni * dzeffddz));
    //      dxnedtt = dxnidtt*zeff+dxnidt*dzeffdt+dxnidt*dzeffdt+xni*dzeffdtt
    dxnedtt = ((((dxnidtt * zeff) + (dxnidt * dzeffdt)) + (dxnidt * dzeffdt)) + (xni * dzeffdtt));
    //      dxnedta = dxnidta*zeff+dxnidt*dzeffda+dxnida*dzeffdt+xni*dzeffdta
    dxnedta = ((((dxnidta * zeff) + (dxnidt * dzeffda)) + (dxnida * dzeffdt)) + (xni * dzeffdta));
    //      dxnedtz = dxnidtz*zeff+dxnidt*dzeffdz+dxnidz*dzeffdt+xni*dzeffdtz
    dxnedtz = ((((dxnidtz * zeff) + (dxnidt * dzeffdz)) + (dxnidz * dzeffdt)) + (xni * dzeffdtz));
    //      dxnedaa = dxnidaa*zeff+dxnida*dzeffda+dxnida*dzeffda+xni*dzeffdaa
    dxnedaa = ((((dxnidaa * zeff) + (dxnida * dzeffda)) + (dxnida * dzeffda)) + (xni * dzeffdaa));
    //      dxnedaz = dxnidaz*zeff+dxnida*dzeffdz+dxnidz*dzeffda+xni*dzeffdaz
    dxnedaz = ((((dxnidaz * zeff) + (dxnida * dzeffdz)) + (dxnidz * dzeffda)) + (xni * dzeffdaz));
    //      dxnedzz = dxnidzz*zeff+dxnidz*dzeffdz+dxnidz*dzeffdz+xni*dzeffdzz
    dxnedzz = ((((dxnidzz * zeff) + (dxnidz * dzeffdz)) + (dxnidz * dzeffdz)) + (xni * dzeffdzz));
    //
    //
    //
    //! first derivatives of the fermi integral electron number densities
    //      dxneferdd = dxnefer_deta * detadd
    dxneferdd = (dxnefer_deta * detadd);
    //      dxneferdt = dxnefer_deta * detadt + dxnefer_dbeta * dbetadt
    dxneferdt = ((dxnefer_deta * detadt) + (dxnefer_dbeta * dbetadt));
    //      dxneferda = dxnefer_deta * detada
    dxneferda = (dxnefer_deta * detada);
    //      dxneferdz = dxnefer_deta * detadz
    dxneferdz = (dxnefer_deta * detadz);
    //
    //! second derivatives
    //      dxneferddd = dxnefer_deta2*detadd*detadd + dxnefer_deta*detaddd
    dxneferddd = (((dxnefer_deta2 * detadd) * detadd) + (dxnefer_deta * detaddd));
    //      dxneferddt = dxnefer_deta2*detadt*detadd &
    //                   + dxnefer_deta_dbeta * dbetadt * detadd &
    //                   + dxnefer_deta *detaddt
    dxneferddt = ((((dxnefer_deta2 * detadt) * detadd) + ((dxnefer_deta_dbeta * dbetadt) * detadd)) +
                  (dxnefer_deta * detaddt));
    //      dxneferdda = dxnefer_deta2*detada*detadd + dxnefer_deta*detadda
    dxneferdda = (((dxnefer_deta2 * detada) * detadd) + (dxnefer_deta * detadda));
    //      dxneferddz = dxnefer_deta2*detadz*detadd + dxnefer_deta*detaddz
    dxneferddz = (((dxnefer_deta2 * detadz) * detadd) + (dxnefer_deta * detaddz));
    //      dxneferdtt = dxnefer_deta2*detadt*detadt + dxnefer_deta*detadtt &
    //                   + 2.0d0*dxnefer_deta_dbeta * detadt * dbetadt &
    //                   + dxnefer_dbeta2 * dbetadt * dbetadt
    dxneferdtt = (((((dxnefer_deta2 * detadt) * detadt) + (dxnefer_deta * detadtt)) +
                   (((2.0 * dxnefer_deta_dbeta) * detadt) * dbetadt)) +
                  ((dxnefer_dbeta2 * dbetadt) * dbetadt));
    //      dxneferdta = dxnefer_deta2 * detada * detadt &
    //                   + dxnefer_deta_dbeta * dbetadt * detada &
    //                   + dxnefer_deta * detadta
    dxneferdta = ((((dxnefer_deta2 * detada) * detadt) + ((dxnefer_deta_dbeta * dbetadt) * detada)) +
                  (dxnefer_deta * detadta));
    //      dxneferdtz = dxnefer_deta2 * detadz * detadt &
    //                   + dxnefer_deta_dbeta * dbetadt * detadz &
    //                   + dxnefer_deta * detadtz
    dxneferdtz = ((((dxnefer_deta2 * detadz) * detadt) + ((dxnefer_deta_dbeta * dbetadt) * detadz)) +
                  (dxnefer_deta * detadtz));
    //      dxneferdaa = dxnefer_deta2*detada*detada + dxnefer_deta*detadaa
    dxneferdaa = (((dxnefer_deta2 * detada) * detada) + (dxnefer_deta * detadaa));
    //      dxneferdaz = dxnefer_deta2*detadz*detada + dxnefer_deta*detadaz
    dxneferdaz = (((dxnefer_deta2 * detadz) * detada) + (dxnefer_deta * detadaz));
    //      dxneferdzz = dxnefer_deta2*detadz*detadz + dxnefer_deta*detadzz
    dxneferdzz = (((dxnefer_deta2 * detadz) * detadz) + (dxnefer_deta * detadzz));
    //
    //
    //
    //! first derivatives of the fermi integral positron number densities
    //      dxnpferdd = dxnpfer_deta * detadd
    dxnpferdd = (dxnpfer_deta * detadd);
    //      dxnpferdt = dxnpfer_deta * detadt + dxnpfer_dbeta * dbetadt
    dxnpferdt = ((dxnpfer_deta * detadt) + (dxnpfer_dbeta * dbetadt));
    //      dxnpferda = dxnpfer_deta * detada
    dxnpferda = (dxnpfer_deta * detada);
    //      dxnpferdz = dxnpfer_deta * detadz
    dxnpferdz = (dxnpfer_deta * detadz);
    //
    //
    //! second derivatives
    //      dxnpferddd = dxnpfer_deta2*detadd*detadd + dxnpfer_deta*detaddd
    dxnpferddd = (((dxnpfer_deta2 * detadd) * detadd) + (dxnpfer_deta * detaddd));
    //      dxnpferddt = dxnpfer_deta2*detadt*detadd &
    //                   + dxnpfer_deta_dbeta * dbetadt * detadd &
    //                   + dxnpfer_deta *detaddt
    dxnpferddt = ((((dxnpfer_deta2 * detadt) * detadd) + ((dxnpfer_deta_dbeta * dbetadt) * detadd)) +
                  (dxnpfer_deta * detaddt));
    //      dxnpferdda = dxnpfer_deta2*detada*detadd + dxnpfer_deta*detadda
    dxnpferdda = (((dxnpfer_deta2 * detada) * detadd) + (dxnpfer_deta * detadda));
    //      dxnpferddz = dxnpfer_deta2*detadz*detadd + dxnpfer_deta*detaddz
    dxnpferddz = (((dxnpfer_deta2 * detadz) * detadd) + (dxnpfer_deta * detaddz));
    //      dxnpferdtt = dxnpfer_deta2*detadt*detadt + dxnpfer_deta*detadtt &
    //                   + 2.0d0*dxnpfer_deta_dbeta * detadt * dbetadt &
    //                   + dxnpfer_dbeta2 * dbetadt * dbetadt
    dxnpferdtt = (((((dxnpfer_deta2 * detadt) * detadt) + (dxnpfer_deta * detadtt)) +
                   (((2.0 * dxnpfer_deta_dbeta) * detadt) * dbetadt)) +
                  ((dxnpfer_dbeta2 * dbetadt) * dbetadt));
    //      dxnpferdta = dxnpfer_deta2 * detada * detadt &
    //                   + dxnpfer_deta_dbeta * dbetadt * detada &
    //                   + dxnpfer_deta * detadta
    dxnpferdta = ((((dxnpfer_deta2 * detada) * detadt) + ((dxnpfer_deta_dbeta * dbetadt) * detada)) +
                  (dxnpfer_deta * detadta));
    //      dxnpferdtz = dxnpfer_deta2 * detadz * detadt &
    //                   + dxnpfer_deta_dbeta * dbetadt * detadz &
    //                   + dxnpfer_deta * detadtz
    dxnpferdtz = ((((dxnpfer_deta2 * detadz) * detadt) + ((dxnpfer_deta_dbeta * dbetadt) * detadz)) +
                  (dxnpfer_deta * detadtz));
    //      dxnpferdaa = dxnpfer_deta2*detada*detada + dxnpfer_deta*detadaa
    dxnpferdaa = (((dxnpfer_deta2 * detada) * detada) + (dxnpfer_deta * detadaa));
    //      dxnpferdaz = dxnpfer_deta2*detadz*detada + dxnpfer_deta*detadaz
    dxnpferdaz = (((dxnpfer_deta2 * detadz) * detada) + (dxnpfer_deta * detadaz));
    //      dxnpferdzz = dxnpfer_deta2*detadz*detadz + dxnpfer_deta*detadzz
    dxnpferdzz = (((dxnpfer_deta2 * detadz) * detadz) + (dxnpfer_deta * detadzz));
    //
    //
    //
    //
    //
    //! now get the pressure and energy
    //! for the electrons
    //
    //      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta, &
    //                                       f32eta2, f32beta2, f32etabeta)
    dfermi(1.5, etaele, beta, f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta);
    //      call dfermi(2.5d0, etaele, beta, f52, f52eta, f52beta, &
    //                                       f52eta2, f52beta2, f52etabeta)
    dfermi(2.5, etaele, beta, f52, f52eta, f52beta, f52eta2, f52beta2, f52etabeta);
    //
    //! pressure in erg/cm**3
    //      yy   = pconst * beta52
    yy = (pconst * beta52);
    //      ww   = pconst * 2.5d0 * beta32
    ww = ((pconst * 2.5) * beta32);
    //      dum1 = f32 + 0.5d0*beta*f52
    dum1 = (f32 + ((0.5 * beta) * f52));
    //      dum2 = f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta
    dum2 = ((f32beta + (0.5 * f52)) + ((0.5 * beta) * f52beta));
    //      dum3 = f32eta + 0.5d0*beta*f52eta
    dum3 = (f32eta + ((0.5 * beta) * f52eta));
    //
    //      pele         = yy * dum1
    pele = (yy * dum1);
    //      dpele_deta   = yy * dum3
    dpele_deta = (yy * dum3);
    //      dpele_dbeta  = ww * dum1 + yy * dum2
    dpele_dbeta = ((ww * dum1) + (yy * dum2));
    //      dpele_deta2  = yy * (f32eta2 + 0.5d0*beta*f52eta2)
    dpele_deta2 = (yy * (f32eta2 + ((0.5 * beta) * f52eta2)));
    //      dpele_dbeta2 = pconst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
    //                     + yy * (f32beta2 + f52beta + 0.5d0*beta*f52beta2)
    dpele_dbeta2 = (((((pconst * beta12) * 3.75) * dum1) + ((2.0 * ww) * dum2)) +
                    (yy * ((f32beta2 + f52beta) + ((0.5 * beta) * f52beta2))));
    //      dpele_deta_dbeta = ww*dum3 + yy*(f32etabeta +0.5d0*f52eta &
    //                     + 0.5d0*beta*f52etabeta)
    dpele_deta_dbeta = ((ww * dum3) + (yy * ((f32etabeta + (0.5 * f52eta)) + ((0.5 * beta) * f52etabeta))));
    //
    //
    //! first derivatives of the electron pressure
    //      dpeledd = dpele_deta * detadd
    dpeledd = (dpele_deta * detadd);
    //      dpeledt = dpele_deta * detadt + dpele_dbeta * dbetadt
    dpeledt = ((dpele_deta * detadt) + (dpele_dbeta * dbetadt));
    //      dpeleda = dpele_deta * detada
    dpeleda = (dpele_deta * detada);
    //      dpeledz = dpele_deta * detadz
    dpeledz = (dpele_deta * detadz);
    //
    //! second derivatives
    //      dpeleddd = dpele_deta2*detadd*detadd + dpele_deta*detaddd
    dpeleddd = (((dpele_deta2 * detadd) * detadd) + (dpele_deta * detaddd));
    //      dpeleddt = dpele_deta2*detadt*detadd &
    //                   + dpele_deta_dbeta * dbetadt * detadd &
    //                   + dpele_deta *detaddt
    dpeleddt = ((((dpele_deta2 * detadt) * detadd) + ((dpele_deta_dbeta * dbetadt) * detadd)) +
                (dpele_deta * detaddt));
    //      dpeledda = dpele_deta2*detada*detadd + dpele_deta*detadda
    dpeledda = (((dpele_deta2 * detada) * detadd) + (dpele_deta * detadda));
    //      dpeleddz = dpele_deta2*detadz*detadd + dpele_deta*detaddz
    dpeleddz = (((dpele_deta2 * detadz) * detadd) + (dpele_deta * detaddz));
    //      dpeledtt = dpele_deta2*detadt*detadt + dpele_deta*detadtt &
    //                   + 2.0d0*dpele_deta_dbeta * detadt * dbetadt &
    //                   + dpele_dbeta2 * dbetadt * dbetadt
    dpeledtt = (((((dpele_deta2 * detadt) * detadt) + (dpele_deta * detadtt)) +
                 (((2.0 * dpele_deta_dbeta) * detadt) * dbetadt)) +
                ((dpele_dbeta2 * dbetadt) * dbetadt));
    //      dpeledta = dpele_deta2 * detada * detadt &
    //                   + dpele_deta_dbeta * dbetadt * detada &
    //                   + dpele_deta * detadta
    dpeledta = ((((dpele_deta2 * detada) * detadt) + ((dpele_deta_dbeta * dbetadt) * detada)) +
                (dpele_deta * detadta));
    //      dpeledtz = dpele_deta2 * detadz * detadt &
    //                   + dpele_deta_dbeta * dbetadt * detadz &
    //                   + dpele_deta * detadtz
    dpeledtz = ((((dpele_deta2 * detadz) * detadt) + ((dpele_deta_dbeta * dbetadt) * detadz)) +
                (dpele_deta * detadtz));
    //      dpeledaa = dpele_deta2*detada*detada + dpele_deta*detadaa
    dpeledaa = (((dpele_deta2 * detada) * detada) + (dpele_deta * detadaa));
    //      dpeledaz = dpele_deta2*detadz*detada + dpele_deta*detadaz
    dpeledaz = (((dpele_deta2 * detadz) * detada) + (dpele_deta * detadaz));
    //      dpeledzz = dpele_deta2*detadz*detadz + dpele_deta*detadzz
    dpeledzz = (((dpele_deta2 * detadz) * detadz) + (dpele_deta * detadzz));
    //
    //
    //! energy in erg/cm**3
    //      yy   = econst * beta52
    yy = (econst * beta52);
    //      ww   = econst * 2.5d0 * beta32
    ww = ((econst * 2.5) * beta32);
    //      dum1 = f32 + beta*f52
    dum1 = (f32 + (beta * f52));
    //      dum2 = f32beta + f52 + beta*f52beta
    dum2 = ((f32beta + f52) + (beta * f52beta));
    //      dum3 = f32eta + beta*f52eta
    dum3 = (f32eta + (beta * f52eta));
    //
    //      eele        = yy * dum1
    eele = (yy * dum1);
    //      deele_deta  = yy * dum3
    deele_deta = (yy * dum3);
    //      deele_dbeta = ww * dum1 + yy * dum2
    deele_dbeta = ((ww * dum1) + (yy * dum2));
    //      deele_deta2  = yy * (f32eta2 + beta*f52eta2)
    deele_deta2 = (yy * (f32eta2 + (beta * f52eta2)));
    //      deele_dbeta2 = econst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
    //                      + yy * (f32beta2 + 2.0d0*f52beta + beta*f52beta2)
    deele_dbeta2 = (((((econst * beta12) * 3.75) * dum1) + ((2.0 * ww) * dum2)) +
                    (yy * ((f32beta2 + (2.0 * f52beta)) + (beta * f52beta2))));
    //      deele_deta_dbeta = ww*dum3 &
    //                      + yy*(f32etabeta + f52eta + beta*f52etabeta)
    deele_deta_dbeta = ((ww * dum3) + (yy * ((f32etabeta + f52eta) + (beta * f52etabeta))));
    //
    //
    //! first derivatives of the electron energy
    //      deeledd = deele_deta * detadd
    deeledd = (deele_deta * detadd);
    //      deeledt = deele_deta * detadt + deele_dbeta * dbetadt
    deeledt = ((deele_deta * detadt) + (deele_dbeta * dbetadt));
    //      deeleda = deele_deta * detada
    deeleda = (deele_deta * detada);
    //      deeledz = deele_deta * detadz
    deeledz = (deele_deta * detadz);
    //
    //! second derivatives
    //      deeleddd = deele_deta2*detadd*detadd + deele_deta*detaddd
    deeleddd = (((deele_deta2 * detadd) * detadd) + (deele_deta * detaddd));
    //      deeleddt = deele_deta2*detadt*detadd &
    //                   + deele_deta_dbeta * dbetadt * detadd &
    //                   + deele_deta *detaddt
    deeleddt = ((((deele_deta2 * detadt) * detadd) + ((deele_deta_dbeta * dbetadt) * detadd)) +
                (deele_deta * detaddt));
    //      deeledda = deele_deta2*detada*detadd + deele_deta*detadda
    deeledda = (((deele_deta2 * detada) * detadd) + (deele_deta * detadda));
    //      deeleddz = deele_deta2*detadz*detadd + deele_deta*detaddz
    deeleddz = (((deele_deta2 * detadz) * detadd) + (deele_deta * detaddz));
    //      deeledtt = deele_deta2*detadt*detadt + deele_deta*detadtt &
    //                   + 2.0d0*deele_deta_dbeta * detadt * dbetadt &
    //                   + deele_dbeta2 * dbetadt * dbetadt
    deeledtt = (((((deele_deta2 * detadt) * detadt) + (deele_deta * detadtt)) +
                 (((2.0 * deele_deta_dbeta) * detadt) * dbetadt)) +
                ((deele_dbeta2 * dbetadt) * dbetadt));
    //      deeledta = deele_deta2 * detada * detadt &
    //                   + deele_deta_dbeta * dbetadt * detada &
    //                   + deele_deta * detadta
    deeledta = ((((deele_deta2 * detada) * detadt) + ((deele_deta_dbeta * dbetadt) * detada)) +
                (deele_deta * detadta));
    //      deeledtz = deele_deta2 * detadz * detadt &
    //                   + deele_deta_dbeta * dbetadt * detadz &
    //                   + deele_deta * detadtz
    deeledtz = ((((deele_deta2 * detadz) * detadt) + ((deele_deta_dbeta * dbetadt) * detadz)) +
                (deele_deta * detadtz));
    //      deeledaa = deele_deta2*detada*detada + deele_deta*detadaa
    deeledaa = (((deele_deta2 * detada) * detada) + (deele_deta * detadaa));
    //      deeledaz = deele_deta2*detadz*detada + deele_deta*detadaz
    deeledaz = (((deele_deta2 * detadz) * detada) + (deele_deta * detadaz));
    //      deeledzz = deele_deta2*detadz*detadz + deele_deta*detadzz
    deeledzz = (((deele_deta2 * detadz) * detadz) + (deele_deta * detadzz));
    //
    //
    //
    //
    //
    //! for the positrons
    //      ppos              = 0.0d0
    ppos = 0.0;
    //      dppos_detap       = 0.0d0
    dppos_detap = 0.0;
    //      dppos_dbeta       = 0.0d0
    dppos_dbeta = 0.0;
    //      dppos_detap2      = 0.0d0
    dppos_detap2 = 0.0;
    //      dppos_dbeta2      = 0.0d0
    dppos_dbeta2 = 0.0;
    //      dppos_detap_dbeta = 0.0d0
    dppos_detap_dbeta = 0.0;
    //      epos              = 0.0d0
    epos = 0.0;
    //      depos_detap       = 0.0d0
    depos_detap = 0.0;
    //      depos_dbeta       = 0.0d0
    depos_dbeta = 0.0;
    //      depos_detap2      = 0.0d0
    depos_detap2 = 0.0;
    //      depos_dbeta2      = 0.0d0
    depos_dbeta2 = 0.0;
    //      depos_detap_dbeta = 0.0d0
    depos_detap_dbeta = 0.0;
    //
    //
    //      if (beta .gt. positron_start) then
    if ((beta > positron_start)) {
        //       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta, &
        //                                        f32eta2, f32beta2, f32etabeta)
        dfermi(1.5, etapos, beta, f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta);
        //       call dfermi(2.5d0, etapos, beta, f52, f52eta, f52beta, &
        //                                        f52eta2, f52beta2, f52etabeta)
        dfermi(2.5, etapos, beta, f52, f52eta, f52beta, f52eta2, f52beta2, f52etabeta);
        //
        //! pressure
        //       yy   = pconst * beta52
        yy = (pconst * beta52);
        //       ww   = pconst * 2.5d0 * beta32
        ww = ((pconst * 2.5) * beta32);
        //       dum1 = f32 + 0.5d0*beta*f52
        dum1 = (f32 + ((0.5 * beta) * f52));
        //       dum2 = f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta
        dum2 = ((f32beta + (0.5 * f52)) + ((0.5 * beta) * f52beta));
        //       dum3 = f32eta + 0.5d0*beta*f52eta
        dum3 = (f32eta + ((0.5 * beta) * f52eta));
        //
        //       ppos         = yy * dum1
        ppos = (yy * dum1);
        //       dppos_detap  = yy * dum3
        dppos_detap = (yy * dum3);
        //       dppos_dbeta  = ww * dum1 + yy * dum2
        dppos_dbeta = ((ww * dum1) + (yy * dum2));
        //       dppos_detap2 = yy * (f32eta2 + 0.5d0*beta*f52eta2)
        dppos_detap2 = (yy * (f32eta2 + ((0.5 * beta) * f52eta2)));
        //       dppos_dbeta2 = pconst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
        //                      + yy * (f32beta2 + f52beta + 0.5d0*beta*f52beta2)
        dppos_dbeta2 = (((((pconst * beta12) * 3.75) * dum1) + ((2.0 * ww) * dum2)) +
                        (yy * ((f32beta2 + f52beta) + ((0.5 * beta) * f52beta2))));
        //       dppos_detap_dbeta = ww*dum3 + yy*(f32etabeta +0.5d0*f52eta &
        //                      + 0.5d0*beta*f52etabeta)
        dppos_detap_dbeta =
            ((ww * dum3) + (yy * ((f32etabeta + (0.5 * f52eta)) + ((0.5 * beta) * f52etabeta))));
        //
        //! energy
        //       yy   = econst * beta52
        yy = (econst * beta52);
        //       ww   = econst * 2.5d0 * beta32
        ww = ((econst * 2.5) * beta32);
        //       dum1 = f32 + beta*f52
        dum1 = (f32 + (beta * f52));
        //       dum2 = f32beta + f52 + beta*f52beta
        dum2 = ((f32beta + f52) + (beta * f52beta));
        //       dum3 = f32eta + beta*f52eta
        dum3 = (f32eta + (beta * f52eta));
        //
        //       epos        = yy * dum1
        epos = (yy * dum1);
        //       depos_detap  = yy * dum3
        depos_detap = (yy * dum3);
        //       depos_dbeta = ww * dum1 + yy * dum2
        depos_dbeta = ((ww * dum1) + (yy * dum2));
        //       depos_detap2  = yy * (f32eta2 + beta*f52eta2)
        depos_detap2 = (yy * (f32eta2 + (beta * f52eta2)));
        //       depos_dbeta2 = econst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
        //                       + yy * (f32beta2 + 2.0d0*f52beta + beta*f52beta2)
        depos_dbeta2 = (((((econst * beta12) * 3.75) * dum1) + ((2.0 * ww) * dum2)) +
                        (yy * ((f32beta2 + (2.0 * f52beta)) + (beta * f52beta2))));
        //       depos_detap_dbeta = ww*dum3 &
        //                       + yy*(f32etabeta + f52eta + beta*f52etabeta)
        depos_detap_dbeta = ((ww * dum3) + (yy * ((f32etabeta + f52eta) + (beta * f52etabeta))));
        //      end if
    }
    //
    //
    //
    //! convert the etap derivatives to eta derivatives
    //! all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta
    //
    //       dppos_deta  = dppos_detap * detap_deta
    dppos_deta = (dppos_detap * detap_deta);
    //       dppos_dbeta = dppos_dbeta + dppos_detap * detap_dbeta
    dppos_dbeta = (dppos_dbeta + (dppos_detap * detap_dbeta));
    //       dppos_deta2 = dppos_detap2 * detap_deta**2 &
    //                       + dppos_detap * detap_deta2
    dppos_deta2 = ((dppos_detap2 * (detap_deta * detap_deta)) + (dppos_detap * detap_deta2));
    //       dppos_dbeta2 = dppos_dbeta2 &
    //                       + 2.0d0 * dppos_detap_dbeta * detap_dbeta &
    //                       + dppos_detap2 * detap_dbeta**2 &
    //                       + dppos_detap * detap_dbeta2
    dppos_dbeta2 = (((dppos_dbeta2 + ((2.0 * dppos_detap_dbeta) * detap_dbeta)) +
                     (dppos_detap2 * (detap_dbeta * detap_dbeta))) +
                    (dppos_detap * detap_dbeta2));
    //       dppos_deta_dbeta = dppos_detap2 * detap_dbeta * detap_deta &
    //                       + dppos_detap_dbeta * detap_deta &
    //                       + dppos_detap * detap_deta_dbeta
    dppos_deta_dbeta = ((((dppos_detap2 * detap_dbeta) * detap_deta) + (dppos_detap_dbeta * detap_deta)) +
                        (dppos_detap * detap_deta_dbeta));
    //
    //! first derivatives of the positron pressure
    //      dpposdd     = dppos_deta * detadd
    dpposdd = (dppos_deta * detadd);
    //      dpposdt     = dppos_deta * detadt + dppos_dbeta * dbetadt
    dpposdt = ((dppos_deta * detadt) + (dppos_dbeta * dbetadt));
    //      dpposda     = dppos_deta * detada
    dpposda = (dppos_deta * detada);
    //      dpposdz     = dppos_deta * detadz
    dpposdz = (dppos_deta * detadz);
    //
    //
    //! second derivatives
    //      dpposddd = dppos_deta2*detadd*detadd + dppos_deta*detaddd
    dpposddd = (((dppos_deta2 * detadd) * detadd) + (dppos_deta * detaddd));
    //      dpposddt = dppos_deta2*detadt*detadd &
    //                   + dppos_deta_dbeta * dbetadt * detadd &
    //                   + dppos_deta *detaddt
    dpposddt = ((((dppos_deta2 * detadt) * detadd) + ((dppos_deta_dbeta * dbetadt) * detadd)) +
                (dppos_deta * detaddt));
    //      dpposdda = dppos_deta2*detada*detadd + dppos_deta*detadda
    dpposdda = (((dppos_deta2 * detada) * detadd) + (dppos_deta * detadda));
    //      dpposddz = dppos_deta2*detadz*detadd + dppos_deta*detaddz
    dpposddz = (((dppos_deta2 * detadz) * detadd) + (dppos_deta * detaddz));
    //      dpposdtt = dppos_deta2*detadt*detadt + dppos_deta*detadtt &
    //                   + 2.0d0*dppos_deta_dbeta * detadt * dbetadt &
    //                   + dppos_dbeta2 * dbetadt * dbetadt
    dpposdtt = (((((dppos_deta2 * detadt) * detadt) + (dppos_deta * detadtt)) +
                 (((2.0 * dppos_deta_dbeta) * detadt) * dbetadt)) +
                ((dppos_dbeta2 * dbetadt) * dbetadt));
    //      dpposdta = dppos_deta2 * detada * detadt &
    //                   + dppos_deta_dbeta * dbetadt * detada &
    //                   + dppos_deta * detadta
    dpposdta = ((((dppos_deta2 * detada) * detadt) + ((dppos_deta_dbeta * dbetadt) * detada)) +
                (dppos_deta * detadta));
    //      dpposdtz = dppos_deta2 * detadz * detadt &
    //                   + dppos_deta_dbeta * dbetadt * detadz &
    //                   + dppos_deta * detadtz
    dpposdtz = ((((dppos_deta2 * detadz) * detadt) + ((dppos_deta_dbeta * dbetadt) * detadz)) +
                (dppos_deta * detadtz));
    //      dpposdaa = dppos_deta2*detada*detada + dppos_deta*detadaa
    dpposdaa = (((dppos_deta2 * detada) * detada) + (dppos_deta * detadaa));
    //      dpposdaz = dppos_deta2*detadz*detada + dppos_deta*detadaz
    dpposdaz = (((dppos_deta2 * detadz) * detada) + (dppos_deta * detadaz));
    //      dpposdzz = dppos_deta2*detadz*detadz + dppos_deta*detadzz
    dpposdzz = (((dppos_deta2 * detadz) * detadz) + (dppos_deta * detadzz));
    //
    //
    //
    //! convert the etap derivatives to eta derivatives
    //! all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta
    //
    //       depos_deta  = depos_detap * detap_deta
    depos_deta = (depos_detap * detap_deta);
    //       depos_dbeta = depos_dbeta + depos_detap * detap_dbeta
    depos_dbeta = (depos_dbeta + (depos_detap * detap_dbeta));
    //       depos_deta2 = depos_detap2 * detap_deta**2 &
    //                       + depos_detap * detap_deta2
    depos_deta2 = ((depos_detap2 * (detap_deta * detap_deta)) + (depos_detap * detap_deta2));
    //       depos_dbeta2 = depos_dbeta2 &
    //                       + 2.0d0 * depos_detap_dbeta * detap_dbeta &
    //                       + depos_detap2 * detap_dbeta**2 &
    //                       + depos_detap * detap_dbeta2
    depos_dbeta2 = (((depos_dbeta2 + ((2.0 * depos_detap_dbeta) * detap_dbeta)) +
                     (depos_detap2 * (detap_dbeta * detap_dbeta))) +
                    (depos_detap * detap_dbeta2));
    //       depos_deta_dbeta = depos_detap2 * detap_dbeta * detap_deta &
    //                       + depos_detap_dbeta * detap_deta &
    //                       + depos_detap * detap_deta_dbeta
    depos_deta_dbeta = ((((depos_detap2 * detap_dbeta) * detap_deta) + (depos_detap_dbeta * detap_deta)) +
                        (depos_detap * detap_deta_dbeta));
    //
    //
    //! first derivatives of the positron energy
    //      deposdd     = depos_deta * detadd
    deposdd = (depos_deta * detadd);
    //      deposdt     = depos_deta * detadt + depos_dbeta * dbetadt
    deposdt = ((depos_deta * detadt) + (depos_dbeta * dbetadt));
    //      deposda     = depos_deta * detada
    deposda = (depos_deta * detada);
    //      deposdz     = depos_deta * detadz
    deposdz = (depos_deta * detadz);
    //
    //! second derivatives
    //      deposddd = depos_deta2*detadd*detadd + depos_deta*detaddd
    deposddd = (((depos_deta2 * detadd) * detadd) + (depos_deta * detaddd));
    //      deposddt = depos_deta2*detadt*detadd &
    //                   + depos_deta_dbeta * dbetadt * detadd &
    //                   + depos_deta *detaddt
    deposddt = ((((depos_deta2 * detadt) * detadd) + ((depos_deta_dbeta * dbetadt) * detadd)) +
                (depos_deta * detaddt));
    //      deposdda = depos_deta2*detada*detadd + depos_deta*detadda
    deposdda = (((depos_deta2 * detada) * detadd) + (depos_deta * detadda));
    //      deposddz = depos_deta2*detadz*detadd + depos_deta*detaddz
    deposddz = (((depos_deta2 * detadz) * detadd) + (depos_deta * detaddz));
    //      deposdtt = depos_deta2*detadt*detadt + depos_deta*detadtt &
    //                   + 2.0d0*depos_deta_dbeta * detadt * dbetadt &
    //                   + depos_dbeta2 * dbetadt * dbetadt
    deposdtt = (((((depos_deta2 * detadt) * detadt) + (depos_deta * detadtt)) +
                 (((2.0 * depos_deta_dbeta) * detadt) * dbetadt)) +
                ((depos_dbeta2 * dbetadt) * dbetadt));
    //      deposdta = depos_deta2 * detada * detadt &
    //                   + depos_deta_dbeta * dbetadt * detada &
    //                   + depos_deta * detadta
    deposdta = ((((depos_deta2 * detada) * detadt) + ((depos_deta_dbeta * dbetadt) * detada)) +
                (depos_deta * detadta));
    //      deposdtz = depos_deta2 * detadz * detadt &
    //                   + depos_deta_dbeta * dbetadt * detadz &
    //                   + depos_deta * detadtz
    deposdtz = ((((depos_deta2 * detadz) * detadt) + ((depos_deta_dbeta * dbetadt) * detadz)) +
                (depos_deta * detadtz));
    //      deposdaa = depos_deta2*detada*detada + depos_deta*detadaa
    deposdaa = (((depos_deta2 * detada) * detada) + (depos_deta * detadaa));
    //      deposdaz = depos_deta2*detadz*detada + depos_deta*detadaz
    deposdaz = (((depos_deta2 * detadz) * detada) + (depos_deta * detadaz));
    //      deposdzz = depos_deta2*detadz*detadz + depos_deta*detadzz
    deposdzz = (((depos_deta2 * detadz) * detadz) + (depos_deta * detadzz));
    //
    //
    //
    //
    //! electron+positron pressure and its derivatives
    //! note: at high temperatures and low densities, dpepdd is very small
    //! and can go negative, so limit it to be positive definite
    //
    //      pep     = pele    + ppos
    pep = (pele + ppos);
    //      dpepdd  = max(dpeledd + dpposdd, 1.0d-30)
    dpepdd = std::max((dpeledd + dpposdd), 1e-30);
    //      dpepdt  = dpeledt + dpposdt
    dpepdt = (dpeledt + dpposdt);
    //      dpepda  = dpeleda + dpposda
    dpepda = (dpeleda + dpposda);
    //      dpepdz  = dpeledz + dpposdz
    dpepdz = (dpeledz + dpposdz);
    //      dpepddd = dpeleddd + dpposddd
    dpepddd = (dpeleddd + dpposddd);
    //      dpepddt = dpeleddt + dpposddt
    dpepddt = (dpeleddt + dpposddt);
    //      dpepdda = dpeledda + dpposdda
    dpepdda = (dpeledda + dpposdda);
    //      dpepddz = dpeleddz + dpposddz
    dpepddz = (dpeleddz + dpposddz);
    //      dpepdtt = dpeledtt + dpposdtt
    dpepdtt = (dpeledtt + dpposdtt);
    //      dpepdta = dpeledta + dpposdta
    dpepdta = (dpeledta + dpposdta);
    //      dpepdtz = dpeledtz + dpposdtz
    dpepdtz = (dpeledtz + dpposdtz);
    //      dpepdaa = dpeledaa + dpposdaa
    dpepdaa = (dpeledaa + dpposdaa);
    //      dpepdaz = dpeledaz + dpposdaz
    dpepdaz = (dpeledaz + dpposdaz);
    //      dpepdzz = dpeledzz + dpposdzz
    dpepdzz = (dpeledzz + dpposdzz);
    //
    //
    //! electron+positron thermal energy and its derivatives
    //      eep     = eele    + epos
    eep = (eele + epos);
    //      deepdd  = deeledd + deposdd
    deepdd = (deeledd + deposdd);
    //      deepdt  = deeledt + deposdt
    deepdt = (deeledt + deposdt);
    //      deepda  = deeleda + deposda
    deepda = (deeleda + deposda);
    //      deepdz  = deeledz + deposdz
    deepdz = (deeledz + deposdz);
    //      deepddd = deeleddd + deposddd
    deepddd = (deeleddd + deposddd);
    //      deepddt = deeleddt + deposddt
    deepddt = (deeleddt + deposddt);
    //      deepdda = deeledda + deposdda
    deepdda = (deeledda + deposdda);
    //      deepddz = deeleddz + deposddz
    deepddz = (deeleddz + deposddz);
    //      deepdtt = deeledtt + deposdtt
    deepdtt = (deeledtt + deposdtt);
    //      deepdta = deeledta + deposdta
    deepdta = (deeledta + deposdta);
    //      deepdtz = deeledtz + deposdtz
    deepdtz = (deeledtz + deposdtz);
    //      deepdaa = deeledaa + deposdaa
    deepdaa = (deeledaa + deposdaa);
    //      deepdaz = deeledaz + deposdaz
    deepdaz = (deeledaz + deposdaz);
    //      deepdzz = deeledzz + deposdzz
    deepdzz = (deeledzz + deposdzz);
    //
    //
    //
    //! electron entropy in erg/gr/kelvin and its derivatives
    //      y       = kerg*deni
    y = (kerg * deni);
    //      sele    = ((pele + eele)*kti - etaele*xnefer) * y
    sele = ((((pele + eele) * kti) - (etaele * xnefer)) * y);
    //
    //
    //! first derivatives
    //      dseledd = ((dpeledd + deeledd)*kti - detadd*xnefer &
    //                  - etaele*dxneferdd)*y - sele*deni
    dseledd =
        ((((((dpeledd + deeledd) * kti) - (detadd * xnefer)) - (etaele * dxneferdd)) * y) - (sele * deni));
    //      dseledt = ((dpeledt + deeledt)*kti - (pele + eele)/(kt*temp) &
    //                  - detadt*xnefer - etaele*dxneferdt)*y
    dseledt = ((((((dpeledt + deeledt) * kti) - ((pele + eele) / (kt * temp))) - (detadt * xnefer)) -
                (etaele * dxneferdt)) *
               y);
    //      dseleda = ((dpeleda + deeleda)*kti - detada*xnefer &
    //                   - etaele*dxneferda)*y
    dseleda = (((((dpeleda + deeleda) * kti) - (detada * xnefer)) - (etaele * dxneferda)) * y);
    //      dseledz = ((dpeledz + deeledz)*kti - detadz*xnefer &
    //                   - etaele*dxneferdz)*y
    dseledz = (((((dpeledz + deeledz) * kti) - (detadz * xnefer)) - (etaele * dxneferdz)) * y);
    //
    //! second derivatives
    //      dseleddd = ((dpeleddd + deeleddd)*kti - detaddd*xnefer &
    //                  - 2.0d0*detadd*dxneferdd - etaele*dxneferddd)*y &
    //                  - 2.0d0*((dpeledd + deeledd)*kti - detadd*xnefer &
    //                      - etaele*dxneferdd)*y*deni &
    //                  + 2.0d0*sele*deni**2
    dseleddd =
        ((((((((dpeleddd + deeleddd) * kti) - (detaddd * xnefer)) - ((2.0 * detadd) * dxneferdd)) -
            (etaele * dxneferddd)) *
           y) -
          (((2.0 * ((((dpeledd + deeledd) * kti) - (detadd * xnefer)) - (etaele * dxneferdd))) * y) * deni)) +
         ((2.0 * sele) * (deni * deni)));
    //      dseleddt = ((dpeleddt + deeleddt)*kti &
    //                   - (dpeledd + deeledd)*kti/temp &
    //                   - detaddt*xnefer - detadd*dxneferdt &
    //                   - detadt*dxneferdd - etaele*dxneferddt)*y &
    //                   - dseledt*deni
    dseleddt =
        (((((((((dpeleddt + deeleddt) * kti) - (((dpeledd + deeledd) * kti) / temp)) - (detaddt * xnefer)) -
             (detadd * dxneferdt)) -
            (detadt * dxneferdd)) -
           (etaele * dxneferddt)) *
          y) -
         (dseledt * deni));
    //      dseledda = ((dpeledda + deeledda)*kti &
    //                   - detadda*xnefer - detadd*dxneferda &
    //                   - detada*dxneferdd - etaele*dxneferdda)*y &
    //                   - dseleda*deni
    dseledda = ((((((((dpeledda + deeledda) * kti) - (detadda * xnefer)) - (detadd * dxneferda)) -
                   (detada * dxneferdd)) -
                  (etaele * dxneferdda)) *
                 y) -
                (dseleda * deni));
    //      dseleddz = ((dpeleddz + deeleddz)*kti &
    //                   - detaddz*xnefer - detadd*dxneferdz &
    //                   - detadz*dxneferdd - etaele*dxneferddz)*y &
    //                   - dseledz*deni
    dseleddz = ((((((((dpeleddz + deeleddz) * kti) - (detaddz * xnefer)) - (detadd * dxneferdz)) -
                   (detadz * dxneferdd)) -
                  (etaele * dxneferddz)) *
                 y) -
                (dseledz * deni));
    //      dseledtt = ((dpeledtt + deeledtt)*kti &
    //                 - 2.0d0*(dpeledt + deeledt)*kti/temp &
    //                 + 2.0d0*(pele + eele)*kti/temp**2 &
    //                 - detadtt*xnefer - 2.0d0*detadt*dxneferdt &
    //                 - etaele*dxneferdtt)*y
    dseledtt = ((((((((dpeledtt + deeledtt) * kti) - (((2.0 * (dpeledt + deeledt)) * kti) / temp)) +
                    (((2.0 * (pele + eele)) * kti) / (temp * temp))) -
                   (detadtt * xnefer)) -
                  ((2.0 * detadt) * dxneferdt)) -
                 (etaele * dxneferdtt)) *
                y);
    //      dseledta = ((dpeledta + deeledta)*kti &
    //                  - (dpeleda + deeleda)*kti/temp &
    //                  - detadta*xnefer - detadt*dxneferda &
    //                  - detada*dxneferdt - etaele*dxneferdta)*y
    dseledta =
        ((((((((dpeledta + deeledta) * kti) - (((dpeleda + deeleda) * kti) / temp)) - (detadta * xnefer)) -
            (detadt * dxneferda)) -
           (detada * dxneferdt)) -
          (etaele * dxneferdta)) *
         y);
    //      dseledtz = ((dpeledtz + deeledtz)*kti &
    //                  - (dpeledz + deeledz)*kti/temp &
    //                  - detadtz*xnefer - detadt*dxneferdz &
    //                  - detadz*dxneferdt - etaele*dxneferdtz)*y
    dseledtz =
        ((((((((dpeledtz + deeledtz) * kti) - (((dpeledz + deeledz) * kti) / temp)) - (detadtz * xnefer)) -
            (detadt * dxneferdz)) -
           (detadz * dxneferdt)) -
          (etaele * dxneferdtz)) *
         y);
    //      dseledaa = ((dpeledaa + deeledaa)*kti - detadaa*xnefer &
    //                 - 2.0d0*detada*dxneferda - etaele*dxneferdaa)*y
    dseledaa = ((((((dpeledaa + deeledaa) * kti) - (detadaa * xnefer)) - ((2.0 * detada) * dxneferda)) -
                 (etaele * dxneferdaa)) *
                y);
    //      dseledaz = ((dpeledaz + deeledaz)*kti - detadaz*xnefer &
    //                   - detada*dxneferdz &
    //                   - detadz*dxneferda  - etaele*dxneferdaz)*y
    dseledaz = (((((((dpeledaz + deeledaz) * kti) - (detadaz * xnefer)) - (detada * dxneferdz)) -
                  (detadz * dxneferda)) -
                 (etaele * dxneferdaz)) *
                y);
    //      dseledzz = ((dpeledzz + deeledzz)*kti - detadzz*xnefer &
    //                 - 2.0d0*detadz*dxneferdz - etaele*dxneferdzz)*y
    dseledzz = ((((((dpeledzz + deeledzz) * kti) - (detadzz * xnefer)) - ((2.0 * detadz) * dxneferdz)) -
                 (etaele * dxneferdzz)) *
                y);
    //
    //
    //
    //! positron entropy in erg/gr/kelvin and its derivatives
    //      spos    = ((ppos + epos)/kt - etapos*xnpfer) * y
    spos = ((((ppos + epos) / kt) - (etapos * xnpfer)) * y);
    //
    //! first derivatives
    //      dsposdd = ((dpposdd + deposdd)*kti &
    //                 - detap_deta*detadd*xnpfer &
    //                 - etapos*dxnpferdd)*y - spos*deni
    dsposdd =
        ((((((dpposdd + deposdd) * kti) - ((detap_deta * detadd) * xnpfer)) - (etapos * dxnpferdd)) * y) -
         (spos * deni));
    //      dsposdt = ((dpposdt + deposdt)*kti - (ppos + epos)/(kt*temp) &
    //                 - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer &
    //                 - etapos*dxnpferdt)*y
    dsposdt = ((((((dpposdt + deposdt) * kti) - ((ppos + epos) / (kt * temp))) -
                 (((detap_deta * detadt) + (detap_dbeta * dbetadt)) * xnpfer)) -
                (etapos * dxnpferdt)) *
               y);
    //      dsposda = ((dpposda + deposda)*kti &
    //                 - detap_deta*detada*xnpfer &
    //                 - etapos*dxnpferda)*y
    dsposda = (((((dpposda + deposda) * kti) - ((detap_deta * detada) * xnpfer)) - (etapos * dxnpferda)) * y);
    //      dsposdz = ((dpposdz + deposdz)*kti - detap_deta*detadz*xnpfer &
    //                   - etapos*dxnpferdz)*y
    dsposdz = (((((dpposdz + deposdz) * kti) - ((detap_deta * detadz) * xnpfer)) - (etapos * dxnpferdz)) * y);
    //
    //! second derivatives
    //      dsposddd = ((dpposddd + deposddd)*kti &
    //                  - detap_deta*detaddd*xnpfer &
    //                  - 2.0d0*detap_deta*detadd*dxnpferdd &
    //                  - etapos*dxnpferddd)*y &
    //                  - 2.0d0*((dpposdd + deposdd)*kti &
    //                  - detap_deta*detadd*xnpfer &
    //                  - etapos*dxnpferdd)*y*deni &
    //                  + 2.0d0*spos*deni**2
    dsposddd = ((((((((dpposddd + deposddd) * kti) - ((detap_deta * detaddd) * xnpfer)) -
                    (((2.0 * detap_deta) * detadd) * dxnpferdd)) -
                   (etapos * dxnpferddd)) *
                  y) -
                 (((2.0 * ((((dpposdd + deposdd) * kti) - ((detap_deta * detadd) * xnpfer)) -
                           (etapos * dxnpferdd))) *
                   y) *
                  deni)) +
                ((2.0 * spos) * (deni * deni)));
    //      dsposddt = ((dpposddt + deposddt)*kti &
    //                   - (dpposdd + deposdd)*kti/temp &
    //                   - detap_deta*detaddt*xnpfer &
    //                   - detap_deta*detadd*dxnpferdt &
    //                   - (detap_deta*detadt + detap_dbeta*dbetadt)*dxnpferdd &
    //                   - etapos*dxnpferddt)*y &
    //                   - dsposdt*deni
    dsposddt = (((((((((dpposddt + deposddt) * kti) - (((dpposdd + deposdd) * kti) / temp)) -
                     ((detap_deta * detaddt) * xnpfer)) -
                    ((detap_deta * detadd) * dxnpferdt)) -
                   (((detap_deta * detadt) + (detap_dbeta * dbetadt)) * dxnpferdd)) -
                  (etapos * dxnpferddt)) *
                 y) -
                (dsposdt * deni));
    //      dsposdda = ((dpposdda + deposdda)*kti &
    //                   - detap_deta*detadda*xnpfer &
    //                   - detap_deta*detadd*dxnpferda &
    //                   - detap_deta*detada*dxnpferdd &
    //                   - etapos*dxnpferdda)*y &
    //                   - dsposda*deni
    dsposdda = ((((((((dpposdda + deposdda) * kti) - ((detap_deta * detadda) * xnpfer)) -
                    ((detap_deta * detadd) * dxnpferda)) -
                   ((detap_deta * detada) * dxnpferdd)) -
                  (etapos * dxnpferdda)) *
                 y) -
                (dsposda * deni));
    //      dsposddz = ((dpposddz + deposddz)*kti &
    //                   - detap_deta*detaddz*xnpfer &
    //                   - detap_deta*detadd*dxnpferdz &
    //                   - detap_deta*detadz*dxnpferdd &
    //                   - etapos*dxnpferddz)*y &
    //                   - dsposdz*deni
    dsposddz = ((((((((dpposddz + deposddz) * kti) - ((detap_deta * detaddz) * xnpfer)) -
                    ((detap_deta * detadd) * dxnpferdz)) -
                   ((detap_deta * detadz) * dxnpferdd)) -
                  (etapos * dxnpferddz)) *
                 y) -
                (dsposdz * deni));
    //
    //      dsposdtt = ((dpposdtt + deposdtt)*kti &
    //                 - 2.0d0*(dpposdt + deposdt)*kti/temp &
    //                 + 2.0d0*(ppos + epos)*kti/temp**2 &
    //                 - (detap_deta2*detadt**2 &
    //                    + 2.0d0*detap_deta_dbeta*detadt*dbetadt &
    //                    + detap_deta*detadtt &
    //                    + detap_dbeta2*dbetadt**2)*xnpfer &
    //                 - 2.0d0*(detap_deta*detadt &
    //                            + detap_dbeta*dbetadt)*dxnpferdt &
    //                 - etapos*dxnpferdtt)*y
    dsposdtt = ((((((((dpposdtt + deposdtt) * kti) - (((2.0 * (dpposdt + deposdt)) * kti) / temp)) +
                    (((2.0 * (ppos + epos)) * kti) / (temp * temp))) -
                   (((((detap_deta2 * (detadt * detadt)) + (((2.0 * detap_deta_dbeta) * detadt) * dbetadt)) +
                      (detap_deta * detadtt)) +
                     (detap_dbeta2 * (dbetadt * dbetadt))) *
                    xnpfer)) -
                  ((2.0 * ((detap_deta * detadt) + (detap_dbeta * dbetadt))) * dxnpferdt)) -
                 (etapos * dxnpferdtt)) *
                y);
    //
    //!      dsposdt = ((dpposdt + deposdt)*kti - (ppos + epos)/(kt*temp)
    //!     1           - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer
    //!     2           - etapos*dxnpferdt)*y
    //
    //      dsposdta = ((dpposdta + deposdta)*kti &
    //                  - (dpposda + deposda)*kti/temp &
    //                  - (detap_deta2*detadt*detada &
    //                      + detap_deta_dbeta*dbetadt*detada &
    //                      + detap_deta*detadta)*xnpfer &
    //                  - (detap_deta*detadt &
    //                            + detap_dbeta*dbetadt)*dxnpferda &
    //                  - detap_deta*detada*dxnpferdt &
    //                  - etapos*dxnpferdta)*y
    dsposdta = ((((((((dpposdta + deposdta) * kti) - (((dpposda + deposda) * kti) / temp)) -
                    (((((detap_deta2 * detadt) * detada) + ((detap_deta_dbeta * dbetadt) * detada)) +
                      (detap_deta * detadta)) *
                     xnpfer)) -
                   (((detap_deta * detadt) + (detap_dbeta * dbetadt)) * dxnpferda)) -
                  ((detap_deta * detada) * dxnpferdt)) -
                 (etapos * dxnpferdta)) *
                y);
    //      dsposdtz = ((dpposdtz + deposdtz)*kti &
    //                  - (dpposdz + deposdz)*kti/temp &
    //                  - (detap_deta2*detadt*detadz &
    //                      + detap_deta_dbeta*dbetadt*detadz &
    //                      + detap_deta*detadtz)*xnpfer &
    //                  - (detap_deta*detadt &
    //                            + detap_dbeta*dbetadt)*dxnpferdz &
    //                  - detap_deta*detadz*dxnpferdt &
    //                  - etapos*dxnpferdtz)*y
    dsposdtz = ((((((((dpposdtz + deposdtz) * kti) - (((dpposdz + deposdz) * kti) / temp)) -
                    (((((detap_deta2 * detadt) * detadz) + ((detap_deta_dbeta * dbetadt) * detadz)) +
                      (detap_deta * detadtz)) *
                     xnpfer)) -
                   (((detap_deta * detadt) + (detap_dbeta * dbetadt)) * dxnpferdz)) -
                  ((detap_deta * detadz) * dxnpferdt)) -
                 (etapos * dxnpferdtz)) *
                y);
    //      dsposdaa = ((dpposdaa + deposdaa)*kti &
    //                  - detap_deta*detadaa*xnpfer &
    //                  - 2.0d0*detap_deta*detada*dxnpferda &
    //                  - etapos*dxnpferdaa)*y
    dsposdaa = ((((((dpposdaa + deposdaa) * kti) - ((detap_deta * detadaa) * xnpfer)) -
                  (((2.0 * detap_deta) * detada) * dxnpferda)) -
                 (etapos * dxnpferdaa)) *
                y);
    //      dsposdaz = ((dpposdaz + deposdaz)*kti &
    //                  - detap_deta*detadaz*xnpfer &
    //                  - detap_deta*detada*dxnpferdz &
    //                  - detap_deta*detadz*dxnpferda &
    //                  - etapos*dxnpferdaz)*y
    dsposdaz = (((((((dpposdaz + deposdaz) * kti) - ((detap_deta * detadaz) * xnpfer)) -
                   ((detap_deta * detada) * dxnpferdz)) -
                  ((detap_deta * detadz) * dxnpferda)) -
                 (etapos * dxnpferdaz)) *
                y);
    //      dsposdzz = ((dpposdzz + deposdzz)*kti &
    //                 - detap_deta*detadzz*xnpfer &
    //                 - 2.0d0*detap_deta*detadz*dxnpferdz &
    //                 - etapos*dxnpferdzz)*y
    dsposdzz = ((((((dpposdzz + deposdzz) * kti) - ((detap_deta * detadzz) * xnpfer)) -
                  (((2.0 * detap_deta) * detadz) * dxnpferdz)) -
                 (etapos * dxnpferdzz)) *
                y);
    //
    //
    //! and their sum
    //      sep      = sele + spos
    sep = (sele + spos);
    //      dsepdd   = dseledd + dsposdd
    dsepdd = (dseledd + dsposdd);
    //      dsepdt   = dseledt + dsposdt
    dsepdt = (dseledt + dsposdt);
    //      dsepda   = dseleda + dsposda
    dsepda = (dseleda + dsposda);
    //      dsepdz   = dseledz + dsposdz
    dsepdz = (dseledz + dsposdz);
    //      dsepddd  = dseleddd + dsposddd
    dsepddd = (dseleddd + dsposddd);
    //      dsepddt  = dseleddt + dsposddt
    dsepddt = (dseleddt + dsposddt);
    //      dsepdda  = dseledda + dsposdda
    dsepdda = (dseledda + dsposdda);
    //      dsepddz  = dseleddz + dsposddz
    dsepddz = (dseleddz + dsposddz);
    //      dsepdtt  = dseledtt + dsposdtt
    dsepdtt = (dseledtt + dsposdtt);
    //      dsepdta  = dseledta + dsposdta
    dsepdta = (dseledta + dsposdta);
    //      dsepdtz  = dseledtz + dsposdtz
    dsepdtz = (dseledtz + dsposdtz);
    //      dsepdaa  = dseledaa + dsposdaa
    dsepdaa = (dseledaa + dsposdaa);
    //      dsepdaz  = dseledaz + dsposdaz
    dsepdaz = (dseledaz + dsposdaz);
    //      dsepdzz  = dseledzz + dsposdzz
    dsepdzz = (dseledzz + dsposdzz);
    //
    //
    //
    //! adjust for the rest mass energy of the positrons
    //      y        = 2.0d0 * mecc
    y = (2.0 * mecc);
    //      epos     = epos     + y * xnpfer
    epos = (epos + (y * xnpfer));
    //      deposdd  = deposdd  + y * dxnpferdd
    deposdd = (deposdd + (y * dxnpferdd));
    //      deposdt  = deposdt  + y * dxnpferdt
    deposdt = (deposdt + (y * dxnpferdt));
    //      deposda  = deposda  + y * dxnpferda
    deposda = (deposda + (y * dxnpferda));
    //      deposdz  = deposdz  + y * dxnpferdz
    deposdz = (deposdz + (y * dxnpferdz));
    //      deposddd = deposddd + y * dxnpferddd
    deposddd = (deposddd + (y * dxnpferddd));
    //      deposddt = deposddt + y * dxnpferddt
    deposddt = (deposddt + (y * dxnpferddt));
    //      deposdda = deposdda + y * dxnpferdda
    deposdda = (deposdda + (y * dxnpferdda));
    //      deposddz = deposddz + y * dxnpferddz
    deposddz = (deposddz + (y * dxnpferddz));
    //      deposdtt = deposdtt + y * dxnpferdtt
    deposdtt = (deposdtt + (y * dxnpferdtt));
    //      deposdta = deposdta + y * dxnpferdta
    deposdta = (deposdta + (y * dxnpferdta));
    //      deposdtz = deposdtz + y * dxnpferdtz
    deposdtz = (deposdtz + (y * dxnpferdtz));
    //      deposdaa = deposdaa + y * dxnpferdaa
    deposdaa = (deposdaa + (y * dxnpferdaa));
    //      deposdaz = deposdaz + y * dxnpferdaz
    deposdaz = (deposdaz + (y * dxnpferdaz));
    //      deposdzz = deposdzz + y * dxnpferdzz
    deposdzz = (deposdzz + (y * dxnpferdzz));
    //
    //
    //! and resum
    //      deepdd  = deeledd + deposdd
    deepdd = (deeledd + deposdd);
    //      deepdt  = deeledt + deposdt
    deepdt = (deeledt + deposdt);
    //      deepda  = deeleda + deposda
    deepda = (deeleda + deposda);
    //      deepdz  = deeledz + deposdz
    deepdz = (deeledz + deposdz);
    //      deepddd = deeleddd + deposddd
    deepddd = (deeleddd + deposddd);
    //      deepddt = deeleddt + deposddt
    deepddt = (deeleddt + deposddt);
    //      deepdda = deeledda + deposdda
    deepdda = (deeledda + deposdda);
    //      deepddz = deeleddz + deposddz
    deepddz = (deeleddz + deposddz);
    //      deepdtt = deeledtt + deposdtt
    deepdtt = (deeledtt + deposdtt);
    //      deepdta = deeledta + deposdta
    deepdta = (deeledta + deposdta);
    //      deepdtz = deeledtz + deposdtz
    deepdtz = (deeledtz + deposdtz);
    //      deepdaa = deeledaa + deposdaa
    deepdaa = (deeledaa + deposdaa);
    //      deepdaz = deeledaz + deposdaz
    deepdaz = (deeledaz + deposdaz);
    //      deepdzz = deeledzz + deposdzz
    deepdzz = (deeledzz + deposdzz);
    //
    //
    //
    //! convert the electron-positron thermal energy in erg/cm**3
    //! to a specific thermal energy in erg/gr
    //
    //      eele     = eele*deni
    eele = (eele * deni);
    //      deeledd  = (deeledd - eele)*deni
    deeledd = ((deeledd - eele) * deni);
    //      deeledt  = deeledt*deni
    deeledt = (deeledt * deni);
    //      deeleda  = deeleda*deni
    deeleda = (deeleda * deni);
    //      deeledz  = deeledz*deni
    deeledz = (deeledz * deni);
    //      deeleddd = (deeleddd - 2.0d0*deeledd)*deni
    deeleddd = ((deeleddd - (2.0 * deeledd)) * deni);
    //      deeleddt = (deeleddt - deeledt)*deni
    deeleddt = ((deeleddt - deeledt) * deni);
    //      deeledda = (deeledda - deeleda)*deni
    deeledda = ((deeledda - deeleda) * deni);
    //      deeleddz = (deeleddz - deeledz)*deni
    deeleddz = ((deeleddz - deeledz) * deni);
    //      deeledtt = deeledtt*deni
    deeledtt = (deeledtt * deni);
    //      deeledta = deeledta*deni
    deeledta = (deeledta * deni);
    //      deeledtz = deeledtz*deni
    deeledtz = (deeledtz * deni);
    //      deeledaa = deeledaa*deni
    deeledaa = (deeledaa * deni);
    //      deeledaz = deeledaz*deni
    deeledaz = (deeledaz * deni);
    //      deeledzz = deeledzz*deni
    deeledzz = (deeledzz * deni);
    //
    //      epos     = epos*deni
    epos = (epos * deni);
    //      deposdd  = (deposdd - epos)*deni
    deposdd = ((deposdd - epos) * deni);
    //      deposdt  = deposdt*deni
    deposdt = (deposdt * deni);
    //      deposda  = deposda*deni
    deposda = (deposda * deni);
    //      deposdz  = deposdz*deni
    deposdz = (deposdz * deni);
    //      deposddd = (deposddd - 2.0d0*deposdd)*deni
    deposddd = ((deposddd - (2.0 * deposdd)) * deni);
    //      deposddt = (deposddt - deposdt)*deni
    deposddt = ((deposddt - deposdt) * deni);
    //      deposdda = (deposdda - deposda)*deni
    deposdda = ((deposdda - deposda) * deni);
    //      deposddz = (deposddz - deposdz)*deni
    deposddz = ((deposddz - deposdz) * deni);
    //      deposdtt = deposdtt*deni
    deposdtt = (deposdtt * deni);
    //      deposdta = deposdta*deni
    deposdta = (deposdta * deni);
    //      deposdtz = deposdtz*deni
    deposdtz = (deposdtz * deni);
    //      deposdaa = deposdaa*deni
    deposdaa = (deposdaa * deni);
    //      deposdaz = deposdaz*deni
    deposdaz = (deposdaz * deni);
    //      deposdzz = deposdzz*deni
    deposdzz = (deposdzz * deni);
    //
    //! and resum
    //      deepdd  = deeledd + deposdd
    deepdd = (deeledd + deposdd);
    //      deepdt  = deeledt + deposdt
    deepdt = (deeledt + deposdt);
    //      deepda  = deeleda + deposda
    deepda = (deeleda + deposda);
    //      deepdz  = deeledz + deposdz
    deepdz = (deeledz + deposdz);
    //      deepddd = deeleddd + deposddd
    deepddd = (deeleddd + deposddd);
    //      deepddt = deeleddt + deposddt
    deepddt = (deeleddt + deposddt);
    //      deepdda = deeledda + deposdda
    deepdda = (deeledda + deposdda);
    //      deepddz = deeleddz + deposddz
    deepddz = (deeleddz + deposddz);
    //      deepdtt = deeledtt + deposdtt
    deepdtt = (deeledtt + deposdtt);
    //      deepdta = deeledta + deposdta
    deepdta = (deeledta + deposdta);
    //      deepdtz = deeledtz + deposdtz
    deepdtz = (deeledtz + deposdtz);
    //      deepdaa = deeledaa + deposdaa
    deepdaa = (deeledaa + deposdaa);
    //      deepdaz = deeledaz + deposdaz
    deepdaz = (deeledaz + deposdaz);
    //      deepdzz = deeledzz + deposdzz
    deepdzz = (deeledzz + deposdzz);
    //
    //
    //
    //! and take care of the ionization potential contributions
    //      if (potmult .eq. 0) then
    if ((potmult == 0)) {
        //       eip     = 0.0d0
        eip = 0.0;
        //       deipdd  = 0.0d0
        deipdd = 0.0;
        //       deipdt  = 0.0d0
        deipdt = 0.0;
        //       deipda  = 0.0d0
        deipda = 0.0;
        //       deipdz  = 0.0d0
        deipdz = 0.0;
        //       deipddd = 0.0d0
        deipddd = 0.0;
        //       deipddt = 0.0d0
        deipddt = 0.0;
        //       deipdda = 0.0d0
        deipdda = 0.0;
        //       deipddz = 0.0d0
        deipddz = 0.0;
        //       deipdtt = 0.0d0
        deipdtt = 0.0;
        //       deipdta = 0.0d0
        deipdta = 0.0;
        //       deipdtz = 0.0d0
        deipdtz = 0.0;
        //       deipdaa = 0.0d0
        deipdaa = 0.0;
        //       deipdaz = 0.0d0
        deipdaz = 0.0;
        //       deipdzz = 0.0d0
        deipdzz = 0.0;
        //
        //       sip     = 0.0d0
        sip = 0.0;
        //       dsipdd  = 0.0d0
        dsipdd = 0.0;
        //       dsipdt  = 0.0d0
        dsipdt = 0.0;
        //       dsipda  = 0.0d0
        dsipda = 0.0;
        //       dsipdz  = 0.0d0
        dsipdz = 0.0;
        //       dsipddd = 0.0d0
        dsipddd = 0.0;
        //       dsipddt = 0.0d0
        dsipddt = 0.0;
        //       dsipdda = 0.0d0
        dsipdda = 0.0;
        //       dsipddz = 0.0d0
        dsipddz = 0.0;
        //       dsipdtt = 0.0d0
        dsipdtt = 0.0;
        //       dsipdta = 0.0d0
        dsipdta = 0.0;
        //       dsipdtz = 0.0d0
        dsipdtz = 0.0;
        //       dsipdaa = 0.0d0
        dsipdaa = 0.0;
        //       dsipdaz = 0.0d0
        dsipdaz = 0.0;
        //       dsipdzz = 0.0d0
        dsipdzz = 0.0;
        //
        //      else
    } else {
        //       eip     = chi * xne
        eip = (chi * xne);
        //       deipdd  = chi * dxnedd
        deipdd = (chi * dxnedd);
        //       deipdt  = chi * dxnedt
        deipdt = (chi * dxnedt);
        //       deipda  = chi * dxneda
        deipda = (chi * dxneda);
        //       deipdz  = chi * dxnedz + hion*ev2erg*xne
        deipdz = ((chi * dxnedz) + ((hion * ev2erg) * xne));
        //       deipddd = chi * dxneddd
        deipddd = (chi * dxneddd);
        //       deipddt = chi * dxneddt
        deipddt = (chi * dxneddt);
        //       deipdda = chi * dxnedda
        deipdda = (chi * dxnedda);
        //       deipddz = chi * dxneddz
        deipddz = (chi * dxneddz);
        //       deipdtt = chi * dxnedtt
        deipdtt = (chi * dxnedtt);
        //       deipdta = chi * dxnedta
        deipdta = (chi * dxnedta);
        //       deipdtz = chi * dxnedtz
        deipdtz = (chi * dxnedtz);
        //       deipdaa = chi * dxnedaa
        deipdaa = (chi * dxnedaa);
        //       deipdaz = chi * dxnedaz
        deipdaz = (chi * dxnedaz);
        //       deipdzz = chi * dxnedzz + 2.0d0*hion*ev2erg*dxnedz
        deipdzz = ((chi * dxnedzz) + (((2.0 * hion) * ev2erg) * dxnedz));
        //
        //
        //! the ionization entropy in erg/gr/kelvin and its derivatives
        //       y       = kerg*deni
        y = (kerg * deni);
        //       sip     = eip*kti*y
        sip = ((eip * kti) * y);
        //       dsipdd  = deipdd*kti*y - sip*deni
        dsipdd = (((deipdd * kti) * y) - (sip * deni));
        //       dsipdt  = (deipdt*kti - eip*kti/temp)*y
        dsipdt = (((deipdt * kti) - ((eip * kti) / temp)) * y);
        //       dsipda  = deipda*kti*y
        dsipda = ((deipda * kti) * y);
        //       dsipdz  = deipdz*kti*y
        dsipdz = ((deipdz * kti) * y);
        //       dsipddd = deipddd*kti*y - dsipdd*deni + sip*deni*deni
        dsipddd = ((((deipddd * kti) * y) - (dsipdd * deni)) + ((sip * deni) * deni));
        //       dsipddt = deipddt*kti*y - dsipdt*deni
        dsipddt = (((deipddt * kti) * y) - (dsipdt * deni));
        //       dsipddt = deipdda*kti*y - dsipda*deni
        dsipddt = (((deipdda * kti) * y) - (dsipda * deni));
        //       dsipddt = deipddz*kti*y - dsipdz*deni
        dsipddt = (((deipddz * kti) * y) - (dsipdz * deni));
        //       dsipdtt = (deipdtt*kti - 2.0d0*deipdt*kti/temp &
        //                   + 2.0d0*eip*kti/temp**2)*y
        dsipdtt =
            ((((deipdtt * kti) - (((2.0 * deipdt) * kti) / temp)) + (((2.0 * eip) * kti) / (temp * temp))) *
             y);
        //       dsipdta = (deipdta*kti - deipda*kti/temp)*y
        dsipdta = (((deipdta * kti) - ((deipda * kti) / temp)) * y);
        //       dsipdtz = (deipdtz*kti - deipdz*kti/temp)*y
        dsipdtz = (((deipdtz * kti) - ((deipdz * kti) / temp)) * y);
        //       dsipdaa = deipdaa*kti*y
        dsipdaa = ((deipdaa * kti) * y);
        //       dsipdaz = deipdaz*kti*y
        dsipdaz = ((deipdaz * kti) * y);
        //       dsipdzz = deipdzz*kti*y
        dsipdzz = ((deipdzz * kti) * y);
        //
        //!       sip    = (eip/kt - etaele*xne) * y
        //!       dsipdd = (deipdd/kt
        //!     1            - detadd*xne)*y
        //!     2            - etaele*dxnedd*y
        //!     3            - sip*deni
        //!       dsipdt = (deipdt/kt
        //!     1             - detadt*xne
        //!     2             - etaele*dxnedt
        //!     3             - eip/(kt*temp))*y
        //
        //! convert the ionization energy from erg/cm**3 to  erg/gr
        //       eip    = eip*deni
        eip = (eip * deni);
        //       deipdd = (deipdd - eip)*deni
        deipdd = ((deipdd - eip) * deni);
        //       deipdt = deipdt*deni
        deipdt = (deipdt * deni);
        //       deipda = deipda*deni
        deipda = (deipda * deni);
        //       deipdz = deipdz*deni
        deipdz = (deipdz * deni);
        //       deipddd = (deipddd - 2.0d0*deipdd)*deni
        deipddd = ((deipddd - (2.0 * deipdd)) * deni);
        //       deipddt = (deipddt - deipdt)*deni
        deipddt = ((deipddt - deipdt) * deni);
        //       deipdda = (deipdda - deipda)*deni
        deipdda = ((deipdda - deipda) * deni);
        //       deipddz = (deipddz - deipdz)*deni
        deipddz = ((deipddz - deipdz) * deni);
        //       deipdtt = deipdtt*deni
        deipdtt = (deipdtt * deni);
        //       deipdta = deipdta*deni
        deipdta = (deipdta * deni);
        //       deipdtz = deipdtz*deni
        deipdtz = (deipdtz * deni);
        //       deipdaa = deipdaa*deni
        deipdaa = (deipdaa * deni);
        //       deipdaz = deipdaz*deni
        deipdaz = (deipdaz * deni);
        //       deipdzz = deipdzz*deni
        deipdzz = (deipdzz * deni);
        //
        //! end of ionization energy block
        //      end if
    }
    //
    //      return
    return;
    //      end
    //
    //
}
void etages(double xni, double zbar, double temp, double &eta) {
    [[maybe_unused]] double xne{};
    [[maybe_unused]] double x{};
    [[maybe_unused]] double y{};
    [[maybe_unused]] double z{};
    [[maybe_unused]] double kt{};
    [[maybe_unused]] double beta{};
    [[maybe_unused]] double tmkt{};
    [[maybe_unused]] double xnefac{};
    [[maybe_unused]] double rt2{};
    [[maybe_unused]] double rt3{};
    [[maybe_unused]] double rtpi{};
    [[maybe_unused]] double cpf0{};
    [[maybe_unused]] double cpf1{};
    [[maybe_unused]] double cpf2{};
    [[maybe_unused]] double cpf3{};
    [[maybe_unused]] double twoth{};
    [[maybe_unused]] double fa0{};
    [[maybe_unused]] double forpi{};
    [[maybe_unused]] double mecc{};
    //      subroutine etages(xni,zbar,temp,eta)
    //      include 'implno.dek'
    //      include 'const.dek'
    //
    //! this routine makes a damn good guess for the electron degeneracy
    //! parameter eta.
    //! input is the ion number density xni, average charge zbar,
    //! and temperature temp.
    //! output is a guess at the electron chemical potential eta
    //
    //
    //! declare the pass
    //      double precision  xni,zbar,temp,eta
    //
    //! declare
    //      double precision xne,x,y,z,kt,beta,tmkt,xnefac
    //
    //      double precision rt2,rt3,rtpi,cpf0,cpf1,cpf2,cpf3, &
    //                       twoth,fa0,forpi,mecc
    //      parameter        (rt2     = 1.4142135623730951d0, &
    //                        rt3     = 1.7320508075688772d0, &
    //                        rtpi    = 1.7724538509055159d0, &
    //                        cpf0    = h/(me*clight), &
    //                        cpf1    = 3.0d0/(8.0d0*pi) * cpf0**3, &
    //                        cpf2    = 4.0d0/cpf1, &
    //                        cpf3    = 2.0d0*rt3*rtpi/(rt2*cpf1), &
    //                        twoth   = 2.0d0/3.0d0, &
    //                        fa0     = 64.0d0/(9.0d0*pi), &
    //                        forpi   = 4.0d0 * pi, &
    //                        mecc    = me * clight * clight)
    rt2 = 1.4142135623730951;
    rt3 = 1.7320508075688772;
    rtpi = 1.7724538509055159;
    cpf0 = (h / (me * clight));
    cpf1 = ((3.0 / (8.0 * pi)) * (cpf0 * cpf0 * cpf0));
    cpf2 = (4.0 / cpf1);
    cpf3 = (((2.0 * rt3) * rtpi) / (rt2 * cpf1));
    twoth = (2.0 / 3.0);
    fa0 = (64.0 / (9.0 * pi));
    forpi = (4.0 * pi);
    mecc = ((me * clight) * clight);
    //
    //! notes: rt2=sqrt(2)  rt3=sqrt(3)  rtpi=sqrt(pi)
    //
    //
    //! for the purposes of guessing eta, assume full ionization
    //      xne   = xni * zbar
    xne = (xni * zbar);
    //      kt    = kerg * temp
    kt = (kerg * temp);
    //      beta  = kt/mecc
    beta = (kt / mecc);
    //
    //
    //! number density of ionized electrons (c&g 24.354k) and number density at
    //! turning point (c&g 24.354i). if either of these exceed the number density
    //! as given by a saha equation, then pairs are important. set eta = -1/2.
    //
    //      if (beta .ge. 1.0) then
    if ((beta >= 1.0)) {
        //       x = cpf2 * beta * beta
        x = ((cpf2 * beta) * beta);
        //      else
    } else {
        //       x = cpf3 * beta * (1.0d0 + 0.75d0*beta) * exp(-1.0d0/beta)
        x = (((cpf3 * beta) * (1.0 + (0.75 * beta))) * std::exp(((-1.0) / beta)));
        //      end if
    }
    //      if (x .ge. xne) then
    if ((x >= xne)) {
        //       eta = -0.5d0
        eta = (-0.5);
        //
        //
        //! get the dimensionless number density (c&g 24.313), if it is large apply the
        //! formula (c&g 24.309) to get a possible alfa, if not large do a two term
        //! binomial expansion on (c&g 24.309) to estimate eta.
        //
        //      else
    } else {
        //       z = (xne*cpf1)**twoth
        z = std::pow((xne * cpf1), twoth);
        //       if (z .ge. 1.0e-6) then
        if ((z >= 1e-06)) {
            //        y = (sqrt(z + 1.0d0) - 1.0d0)/beta
            y = ((std::sqrt((z + 1.0)) - 1.0) / beta);
            //       else
        } else {
            //        y = z * (1.0d0 - z * 0.25d0) * 0.5d0/beta
            y = (((z * (1.0 - (z * 0.25))) * 0.5) / beta);
            //       end if
        }
        //
        //
        //! isolate the constant in front of the number density integral. if it is
        //! small enough run the divine approximation backwards with c&g 24.43. then
        //! join it smoothly with the lower limit.
        //
        //       x = log10(xne**0.6d0/temp)
        x = std::log10((std::pow(xne, 0.6) / temp));
        //       if (x .le. 9.5) then
        if ((x <= 9.5)) {
            //        z = ((1.0d0 + fa0*beta)*sqrt(1.0d0 + fa0*beta*0.5) - 1.0d0)/fa0
            z = ((((1.0 + (fa0 * beta)) * std::sqrt((1.0 + ((fa0 * beta) * 0.5)))) - 1.0) / fa0);
            //        tmkt    = 2.0d0 * me/h * kt/h
            tmkt = ((((2.0 * me) / h) * kt) / h);
            //        xnefac  = forpi * tmkt * sqrt(tmkt)
            xnefac = ((forpi * tmkt) * std::sqrt(tmkt));
            //        eta = -log(xnefac*rtpi*(0.5d0+0.75d0*z)/xne)
            eta = (-std::log((((xnefac * rtpi) * (0.5 + (0.75 * z))) / xne)));
            //        if (x .ge. 8.5) eta = eta*(9.5d0-x) + y * (1.0d0 - (9.5d0-x))
            if ((x >= 8.5))
                eta = ((eta * (9.5 - x)) + (y * (1.0 - (9.5 - x))));
            //       else
        } else {
            //        eta = y
            eta = y;
            //       end if
        }
        //      end if
    }
    //
    //      return
    return;
    //      end
    //
    //
    //
    //
    //
}
void dfermi(double dk, double eta, double theta, double &fd, double &fdeta, double &fdtheta, double &fdeta2,
            double &fdtheta2, double &fdetadtheta) {
    [[maybe_unused]] double d{};
    [[maybe_unused]] double sg{};
    [[maybe_unused]] double a1{};
    [[maybe_unused]] double b1{};
    [[maybe_unused]] double c1{};
    [[maybe_unused]] double a2{};
    [[maybe_unused]] double b2{};
    [[maybe_unused]] double c2{};
    [[maybe_unused]] double d2{};
    [[maybe_unused]] double e2{};
    [[maybe_unused]] double a3{};
    [[maybe_unused]] double b3{};
    [[maybe_unused]] double c3{};
    [[maybe_unused]] double d3{};
    [[maybe_unused]] double e3{};
    [[maybe_unused]] double eta1{};
    [[maybe_unused]] double xi{};
    [[maybe_unused]] double xi2{};
    [[maybe_unused]] double x1{};
    [[maybe_unused]] double x2{};
    [[maybe_unused]] double x3{};
    [[maybe_unused]] double s1{};
    [[maybe_unused]] double s2{};
    [[maybe_unused]] double s3{};
    [[maybe_unused]] double s12{};
    double par[4]{};
    double res[5]{};
    double drde[5]{};
    double drdt[5]{};
    double drde2[5]{};
    double drdt2[5]{};
    double drdet[5]{};
    //      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
    //                        fdeta2,fdtheta2,fdetadtheta)
    //      include 'implno.dek'
    //
    //! this routine computes the fermi-dirac integrals f(dk,eta,theta) of
    //! index dk, with degeneracy parameter eta and relativity parameter theta.
    //
    //! input is dk the double precision index of the fermi-dirac function,
    //! eta the degeneracy parameter, and theta the relativity parameter.
    //
    //!. output is fd is computed by applying three 10-point gauss-legendre
    //! and one 10-point gauss-laguerre rules over four appropriate subintervals.
    //! fdeta is the derivative of f with respect to eta,
    //! fdthetha is the derivative of f with theta.
    //
    //! reference: j.m. aparicio, apjs 117, 632 1998
    //
    //
    //! declare
    //      external         fdfunc1,fdfunc2
    //      double precision dk,eta,theta,fd,fdeta,fdtheta, &
    //                       fdeta2,fdtheta2,fdetadtheta, &
    //                       d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3, &
    //                       eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3), &
    //                       res(4),drde(4),drdt(4),drde2(4),drdt2(4),drdet(4)
    //
    //
    //!   parameters defining the location of the breakpoints for the
    //!   subintervals of integration:
    //      data d   / 3.3609d0 /
    d = 3.3609;
    //      data sg  / 9.1186d-2 /
    sg = 0.091186;
    //      data a1  / 6.7774d0 /
    a1 = 6.7774;
    //      data b1  / 1.1418d0 /
    b1 = 1.1418;
    //      data c1  / 2.9826d0 /
    c1 = 2.9826;
    //      data a2  / 3.7601d0 /
    a2 = 3.7601;
    //      data b2  / 9.3719d-2 /
    b2 = 0.093719;
    //      data c2  / 2.1063d-2 /
    c2 = 0.021063;
    //      data d2  / 3.1084d1 /
    d2 = 31.084;
    //      data e2  / 1.0056d0 /
    e2 = 1.0056;
    //      data a3  / 7.5669d0 /
    a3 = 7.5669;
    //      data b3  / 1.1695d0 /
    b3 = 1.1695;
    //      data c3  / 7.5416d-1 /
    c3 = 0.75416;
    //      data d3  / 6.6558d0 /
    d3 = 6.6558;
    //      data e3  /-1.2819d-1 /
    e3 = (-0.12819);
    //
    //
    //!   integrand parameters:
    //      par(1)=dk
    par[1] = dk;
    //      par(2)=eta
    par[2] = eta;
    //      par(3)=theta
    par[3] = theta;
    //
    //
    //!   definition of xi:
    //      eta1=sg*(eta-d)
    eta1 = (sg * (eta - d));
    //      if (eta1.le.5.d1) then
    if ((eta1 <= 50.0)) {
        //        xi=log(1.d0+exp(eta1))/sg
        xi = (std::log((1.0 + std::exp(eta1))) / sg);
        //      else
    } else {
        //        xi=eta-d
        xi = (eta - d);
        //      endif
    }
    //      xi2=xi*xi
    xi2 = (xi * xi);
    //
    //!   definition of the x_i:
    //      x1=(a1  +b1*xi+c1*   xi2) &
    //        /(1.d0+c1*xi)
    x1 = (((a1 + (b1 * xi)) + (c1 * xi2)) / (1.0 + (c1 * xi)));
    //      x2=(a2  +b2*xi+c2*d2*xi2) &
    //        /(1.d0+e2*xi+c2*   xi2)
    x2 = (((a2 + (b2 * xi)) + ((c2 * d2) * xi2)) / ((1.0 + (e2 * xi)) + (c2 * xi2)));
    //      x3=(a3  +b3*xi+c3*d3*xi2) &
    //        /(1.d0+e3*xi+c3*   xi2)
    x3 = (((a3 + (b3 * xi)) + ((c3 * d3) * xi2)) / ((1.0 + (e3 * xi)) + (c3 * xi2)));
    //
    //!   breakpoints:
    //      s1=x1-x2
    s1 = (x1 - x2);
    //      s2=x1
    s2 = x1;
    //      s3=x1+x3
    s3 = (x1 + x3);
    //      s12=sqrt(s1)
    s12 = std::sqrt(s1);
    //
    //!   quadrature integrations:
    //
    //! 14 significant figure accuracy
    //      call dqleg020(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
    //                    drde2(1), drdt2(1), drdet(1), par,3)
    dqleg020(fdfunc2, 0.0, s12, res[1], drde[1], drdt[1], drde2[1], drdt2[1], drdet[1], par, 3);
    //      call dqleg020(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
    //                    drde2(2), drdt2(2), drdet(2), par,3)
    dqleg020(fdfunc1, s1, s2, res[2], drde[2], drdt[2], drde2[2], drdt2[2], drdet[2], par, 3);
    //      call dqleg020(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
    //                    drde2(3), drdt2(3), drdet(3), par,3)
    dqleg020(fdfunc1, s2, s3, res[3], drde[3], drdt[3], drde2[3], drdt2[3], drdet[3], par, 3);
    //      call dqlag020(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
    //                    drde2(4), drdt2(4), drdet(4), par,3)
    dqlag020(fdfunc1, s3, 1.0, res[4], drde[4], drdt[4], drde2[4], drdt2[4], drdet[4], par, 3);
    //
    //! sum the contributions
    //      fd          = res(1)   + res(2)   + res(3)   + res(4)
    fd = (((res[1] + res[2]) + res[3]) + res[4]);
    //      fdeta       = drde(1)  + drde(2)  + drde(3)  + drde(4)
    fdeta = (((drde[1] + drde[2]) + drde[3]) + drde[4]);
    //      fdtheta     = drdt(1)  + drdt(2)  + drdt(3)  + drdt(4)
    fdtheta = (((drdt[1] + drdt[2]) + drdt[3]) + drdt[4]);
    //      fdeta2      = drde2(1) + drde2(2) + drde2(3) + drde2(4)
    fdeta2 = (((drde2[1] + drde2[2]) + drde2[3]) + drde2[4]);
    //      fdtheta2    = drdt2(1) + drdt2(2) + drdt2(3) + drdt2(4)
    fdtheta2 = (((drdt2[1] + drdt2[2]) + drdt2[3]) + drdt2[4]);
    //      fdetadtheta = drdet(1) + drdet(2) + drdet(3) + drdet(4)
    fdetadtheta = (((drdet[1] + drdet[2]) + drdet[3]) + drdet[4]);
    //      return
    return;
    //      end
    //
    //
    //
}
void fdfunc1(double x, const double *par, [[maybe_unused]] int n, double &fd, double &fdeta, double &fdtheta, double &fdeta2,
             double &fdtheta2, double &fdetadtheta) {
    [[maybe_unused]] double dk{};
    [[maybe_unused]] double eta{};
    [[maybe_unused]] double theta{};
    [[maybe_unused]] double factor{};
    [[maybe_unused]] double xst{};
    [[maybe_unused]] double dxst{};
    [[maybe_unused]] double denom{};
    [[maybe_unused]] double denomi{};
    [[maybe_unused]] double denom2{};
    [[maybe_unused]] double xdk{};
    //      subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta, &
    //                         fdeta2,fdtheta2,fdetadtheta)
    //      include 'implno.dek'
    //
    //! forms the fermi-dirac integrand and its first and second
    //! derivatives with eta and theta.
    //
    //! input:
    //! x is the integration variable
    //! par(1) is the double precision index
    //! par(2) is the degeneravy parameter
    //! par(3) is the relativity parameter
    //
    //! output:
    //! fd is the integrand
    //! fdeta is the first derivative with eta
    //! fdeta2 is the second derivative with eta
    //! fdtheta is the first derivative with theta
    //! fdtheta2 is the second derivative with theta
    //! fdetadtheta is the mixed second derivative
    //
    //! declare the pass
    //      integer          n
    //      double precision x,par(n),fd, &
    //                       fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta
    //
    //! local variables
    //      double precision dk,eta,theta, &
    //                       factor,xst,dxst,denom,denomi,denom2,xdk
    //
    //
    //! initialize
    //      dk    = par(1)
    dk = par[1];
    //      eta   = par(2)
    eta = par[2];
    //      theta = par(3)
    theta = par[3];
    //      xdk   = x**dk
    xdk = std::pow(x, dk);
    //      xst   = 1.0d0 + 0.5d0*x*theta
    xst = (1.0 + ((0.5 * x) * theta));
    //      dxst  = sqrt(xst)
    dxst = std::sqrt(xst);
    //
    //!   avoid overflow in the exponentials at large x
    //      if ((x-eta) .lt. 1.0d2) then
    if (((x - eta) < 100.0)) {
        //       factor      = exp(x-eta)
        factor = std::exp((x - eta));
        //       denom       = factor + 1.0d0
        denom = (factor + 1.0);
        //       denomi      = 1.0d0/denom
        denomi = (1.0 / denom);
        //       fd          = xdk * dxst * denomi
        fd = ((xdk * dxst) * denomi);
        //       fdeta       = fd * factor * denomi
        fdeta = ((fd * factor) * denomi);
        //       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
        fdeta2 = ((((2.0 * factor) * denomi) - 1.0) * fdeta);
        //       denom2      = 1.0d0/(4.0d0 * xst)
        denom2 = (1.0 / (4.0 * xst));
        //       fdtheta     = fd * x * denom2
        fdtheta = ((fd * x) * denom2);
        //       fdtheta2    = -fdtheta * x * denom2
        fdtheta2 = (((-fdtheta) * x) * denom2);
        //       fdetadtheta = fdtheta * factor * denomi
        fdetadtheta = ((fdtheta * factor) * denomi);
        //      else
    } else {
        //       factor      = exp(eta-x)
        factor = std::exp((eta - x));
        //       fd          = xdk * dxst * factor
        fd = ((xdk * dxst) * factor);
        //       fdeta       = fd
        fdeta = fd;
        //       fdeta2      = fd
        fdeta2 = fd;
        //       denom2      = 1.0d0/(4.0d0 * xst)
        denom2 = (1.0 / (4.0 * xst));
        //       fdtheta     = fd * x * denom2
        fdtheta = ((fd * x) * denom2);
        //       fdtheta2    = -fdtheta * x * denom2
        fdtheta2 = (((-fdtheta) * x) * denom2);
        //       fdetadtheta = fdtheta
        fdetadtheta = fdtheta;
        //      endif
    }
    //
    //      return
    return;
    //      end
    //
    //
    //
}
void fdfunc2(double x, const double *par, [[maybe_unused]] int n, double &fd, double &fdeta, double &fdtheta, double &fdeta2,
             double &fdtheta2, double &fdetadtheta) {
    [[maybe_unused]] double dk{};
    [[maybe_unused]] double eta{};
    [[maybe_unused]] double theta{};
    [[maybe_unused]] double factor{};
    [[maybe_unused]] double xst{};
    [[maybe_unused]] double dxst{};
    [[maybe_unused]] double denom{};
    [[maybe_unused]] double denomi{};
    [[maybe_unused]] double denom2{};
    [[maybe_unused]] double xdk{};
    [[maybe_unused]] double xsq{};
    //      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta, &
    //                         fdeta2,fdtheta2,fdetadtheta)
    //      include 'implno.dek'
    //
    //! forms the fermi-dirac integrand and its first and second
    //! derivatives with eta and theta when the variable change z**2=x
    //! has been made.
    //
    //! input:
    //! x is the integration variable
    //! par(1) is the double precision index
    //! par(2) is the degeneravy parameter
    //! par(3) is the relativity parameter
    //
    //! output:
    //! fd is the integrand
    //! fdeta is the first derivative with eta
    //! fdeta2 is the second derivative with eta
    //! fdtheta is the first derivative with theta
    //! fdtheta2 is the second derivative with theta
    //! fdetadtheta is the mixed second derivative
    //
    //! declare the pass
    //      integer          n
    //      double precision x,par(n),fd, &
    //                       fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta
    //
    //! local variables
    //      double precision dk,eta,theta, &
    //                       factor,xst,dxst,denom,denomi,denom2,xdk,xsq
    //
    //
    //      dk    = par(1)
    dk = par[1];
    //      eta   = par(2)
    eta = par[2];
    //      theta = par(3)
    theta = par[3];
    //      xsq   = x * x
    xsq = (x * x);
    //      xdk   = x**(2.0d0 * dk + 1.0d0)
    xdk = std::pow(x, ((2.0 * dk) + 1.0));
    //      xst   = 1.0d0 + 0.5d0 * xsq * theta
    xst = (1.0 + ((0.5 * xsq) * theta));
    //      dxst  = sqrt(xst)
    dxst = std::sqrt(xst);
    //
    //!   avoid an overflow in the denominator at large x:
    //      if ((xsq-eta) .lt. 1.d2) then
    if (((xsq - eta) < 100.0)) {
        //       factor      = exp(xsq - eta)
        factor = std::exp((xsq - eta));
        //       denom       = factor + 1.0d0
        denom = (factor + 1.0);
        //       denomi      = 1.0d0/denom
        denomi = (1.0 / denom);
        //       fd          = 2.0d0 * xdk * dxst * denomi
        fd = (((2.0 * xdk) * dxst) * denomi);
        //       fdeta       = fd * factor * denomi
        fdeta = ((fd * factor) * denomi);
        //       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
        fdeta2 = ((((2.0 * factor) * denomi) - 1.0) * fdeta);
        //       denom2      = 1.0d0/(4.0d0 * xst)
        denom2 = (1.0 / (4.0 * xst));
        //       fdtheta     = fd * xsq * denom2
        fdtheta = ((fd * xsq) * denom2);
        //       fdtheta2    = -fdtheta * xsq * denom2
        fdtheta2 = (((-fdtheta) * xsq) * denom2);
        //       fdetadtheta = fdtheta * factor * denomi
        fdetadtheta = ((fdtheta * factor) * denomi);
        //      else
    } else {
        //       factor      = exp(eta - xsq)
        factor = std::exp((eta - xsq));
        //       fd          = 2.0d0 * xdk * dxst * factor
        fd = (((2.0 * xdk) * dxst) * factor);
        //       fdeta       = fd
        fdeta = fd;
        //       fdeta2      = fd
        fdeta2 = fd;
        //       denom2      = 1.0d0/(4.0d0 * xst)
        denom2 = (1.0 / (4.0 * xst));
        //       fdtheta     = fd * xsq * denom2
        fdtheta = ((fd * xsq) * denom2);
        //       fdtheta2    = -fdtheta * xsq * denom2
        fdtheta2 = (((-fdtheta) * xsq) * denom2);
        //       fdetadtheta = fdtheta
        fdetadtheta = fdtheta;
        //      endif
    }
    //
    //      return
    return;
    //      end
    //
    //
    //
}
void dqleg020(Integrand f, double a, double b, double &result, double &drdeta, double &drdtheta,
              double &drdeta2, double &drdtheta2, double &drdetadtheta, const double *par, int n) {
    [[maybe_unused]] int j{};
    [[maybe_unused]] double absc1{};
    [[maybe_unused]] double absc2{};
    [[maybe_unused]] double center{};
    [[maybe_unused]] double hlfrun{};
    double wg[11]{};
    double xg[11]{};
    double fval[3]{};
    double deta[3]{};
    double dtheta[3]{};
    double deta2[3]{};
    double dtheta2[3]{};
    double detadtheta[3]{};
    //      subroutine dqleg020(f,a,b,result,drdeta,drdtheta, &
    //                          drdeta2,drdtheta2,drdetadtheta,par,n)
    //      include 'implno.dek'
    //!
    //! 20 point gauss-legendre rule for the fermi-dirac function and
    //! its derivatives with respect to eta and theta.
    //! on input f is the name of the subroutine containing the integrand,
    //! a is the lower end point of the interval, b is the higher end point,
    //! par is an array of constant parameters to be passed to subroutine f,
    //! and n is the length of the par array. on output result is the
    //! approximation from applying the 20-point gauss-legendre rule,
    //! drdeta is the derivative with respect to eta, and drdtheta is the
    //! derivative with respect to theta.
    //!
    //! note: since the number of nodes is even, zero is not an abscissa.
    //!
    //! declare
    //      external         f
    //      integer          j,n
    //      double precision a,b,result,drdeta,drdtheta, &
    //                       drdeta2,drdtheta2,drdetadtheta,par(n), &
    //                       absc1,absc2,center,hlfrun,wg(10),xg(10), &
    //                       fval(2),deta(2),dtheta(2), &
    //                       deta2(2),dtheta2(2),detadtheta(2)
    //
    //
    //! the abscissae and weights are given for the interval (-1,1).
    //! xg     - abscissae of the 20-point gauss-legendre rule
    //!          for half of the usual run (-1,1), i.e.
    //!          the positive nodes of the 20-point rule
    //! wg     - weights of the 20-point gauss rule.
    //!
    //! abscissae and weights were evaluated with 100 decimal digit arithmetic.
    //
    //
    //      data xg (  1) /   7.65265211334973337546404093988382110d-2 /
    xg[1] = 0.07652652113349734;
    //      data xg (  2) /   2.27785851141645078080496195368574624d-1 /
    xg[2] = 0.22778585114164507;
    //      data xg (  3) /   3.73706088715419560672548177024927237d-1 /
    xg[3] = 0.37370608871541955;
    //      data xg (  4) /   5.10867001950827098004364050955250998d-1 /
    xg[4] = 0.5108670019508271;
    //      data xg (  5) /   6.36053680726515025452836696226285936d-1 /
    xg[5] = 0.636053680726515;
    //      data xg (  6) /   7.46331906460150792614305070355641590d-1 /
    xg[6] = 0.7463319064601508;
    //      data xg (  7) /   8.39116971822218823394529061701520685d-1 /
    xg[7] = 0.8391169718222188;
    //      data xg (  8) /   9.12234428251325905867752441203298113d-1 /
    xg[8] = 0.912234428251326;
    //      data xg (  9) /   9.63971927277913791267666131197277221d-1 /
    xg[9] = 0.9639719272779138;
    //      data xg ( 10) /   9.93128599185094924786122388471320278d-1 /
    xg[10] = 0.9931285991850949;
    //
    //      data wg (  1) /   1.52753387130725850698084331955097593d-1 /
    wg[1] = 0.15275338713072584;
    //      data wg (  2) /   1.49172986472603746787828737001969436d-1 /
    wg[2] = 0.14917298647260374;
    //      data wg (  3) /   1.42096109318382051329298325067164933d-1 /
    wg[3] = 0.14209610931838204;
    //      data wg (  4) /   1.31688638449176626898494499748163134d-1 /
    wg[4] = 0.13168863844917664;
    //      data wg (  5) /   1.18194531961518417312377377711382287d-1 /
    wg[5] = 0.11819453196151841;
    //      data wg (  6) /   1.01930119817240435036750135480349876d-1 /
    wg[6] = 0.10193011981724044;
    //      data wg (  7) /   8.32767415767047487247581432220462061d-2 /
    wg[7] = 0.08327674157670475;
    //      data wg (  8) /   6.26720483341090635695065351870416063d-2 /
    wg[8] = 0.06267204833410907;
    //      data wg (  9) /   4.06014298003869413310399522749321098d-2 /
    wg[9] = 0.04060142980038694;
    //      data wg ( 10) /   1.76140071391521183118619623518528163d-2 /
    wg[10] = 0.017614007139152118;
    //
    //
    //!           list of major variables
    //!           -----------------------
    //!
    //!           absc   - abscissa
    //!           fval*  - function value
    //!           result - result of the 20-point gauss formula
    //
    //
    //      center       = 0.5d0 * (a+b)
    center = (0.5 * (a + b));
    //      hlfrun       = 0.5d0 * (b-a)
    hlfrun = (0.5 * (b - a));
    //      result       = 0.0d0
    result = 0.0;
    //      drdeta       = 0.0d0
    drdeta = 0.0;
    //      drdtheta     = 0.0d0
    drdtheta = 0.0;
    //      drdeta2      = 0.0d0
    drdeta2 = 0.0;
    //      drdtheta2    = 0.0d0
    drdtheta2 = 0.0;
    //      drdetadtheta = 0.0d0
    drdetadtheta = 0.0;
    //      do j=1,10
    for (j = 1; j <= 10; ++j) {
        //        absc1 = center + hlfrun*xg(j)
        absc1 = (center + (hlfrun * xg[j]));
        //        absc2 = center - hlfrun*xg(j)
        absc2 = (center - (hlfrun * xg[j]));
        //
        //        call f(absc1, par, n, fval(1), deta(1), dtheta(1), &
        //               deta2(1), dtheta2(1), detadtheta(1))
        f(absc1, par, n, fval[1], deta[1], dtheta[1], deta2[1], dtheta2[1], detadtheta[1]);
        //        call f(absc2, par, n, fval(2), deta(2), dtheta(2), &
        //               deta2(2), dtheta2(2), detadtheta(2))
        f(absc2, par, n, fval[2], deta[2], dtheta[2], deta2[2], dtheta2[2], detadtheta[2]);
        //
        //        result       = result    + (fval(1)    + fval(2))*wg(j)
        result = (result + ((fval[1] + fval[2]) * wg[j]));
        //        drdeta       = drdeta    + (deta(1)    + deta(2))*wg(j)
        drdeta = (drdeta + ((deta[1] + deta[2]) * wg[j]));
        //        drdtheta     = drdtheta  + (dtheta(1)  + dtheta(2))*wg(j)
        drdtheta = (drdtheta + ((dtheta[1] + dtheta[2]) * wg[j]));
        //        drdeta2      = drdeta2   + (deta2(1)   + deta2(2))*wg(j)
        drdeta2 = (drdeta2 + ((deta2[1] + deta2[2]) * wg[j]));
        //        drdtheta2    = drdtheta2 + (dtheta2(1) + dtheta2(2))*wg(j)
        drdtheta2 = (drdtheta2 + ((dtheta2[1] + dtheta2[2]) * wg[j]));
        //        drdetadtheta = drdetadtheta+(detadtheta(1)+detadtheta(2))*wg(j)
        drdetadtheta = (drdetadtheta + ((detadtheta[1] + detadtheta[2]) * wg[j]));
        //      enddo
    }
    //
    //      result       = result * hlfrun
    result = (result * hlfrun);
    //      drdeta       = drdeta * hlfrun
    drdeta = (drdeta * hlfrun);
    //      drdtheta     = drdtheta * hlfrun
    drdtheta = (drdtheta * hlfrun);
    //      drdeta2      = drdeta2 * hlfrun
    drdeta2 = (drdeta2 * hlfrun);
    //      drdtheta2    = drdtheta2 * hlfrun
    drdtheta2 = (drdtheta2 * hlfrun);
    //      drdetadtheta = drdetadtheta * hlfrun
    drdetadtheta = (drdetadtheta * hlfrun);
    //      return
    return;
    //      end
    //
    //
    //
    //
    //
    //
}
void dqlag020(Integrand f, double a, double b, double &result, double &drdeta, double &drdtheta,
              double &drdeta2, double &drdtheta2, double &drdetadtheta, const double *par, int n) {
    [[maybe_unused]] int j{};
    [[maybe_unused]] double absc{};
    double wg[21]{};
    double xg[21]{};
    [[maybe_unused]] double fval{};
    [[maybe_unused]] double deta{};
    [[maybe_unused]] double dtheta{};
    [[maybe_unused]] double deta2{};
    [[maybe_unused]] double dtheta2{};
    [[maybe_unused]] double detadtheta{};
    //      subroutine dqlag020(f,a,b,result,drdeta,drdtheta, &
    //                          drdeta2,drdtheta2,drdetadtheta,par,n)
    //      include 'implno.dek'
    //!
    //! 20 point gauss-laguerre rule for the fermi-dirac function.
    //! on input f is the external function defining the integrand
    //! f(x)=g(x)*w(x), where w(x) is the gaussian weight
    //! w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
    //! a is the lower end point of the interval, b is the higher end point,
    //! par is an array of constant parameters to be passed to the function f,
    //! and n is the length of the par array. on output result is the
    //! approximation from applying the 20-point gauss-laguerre rule.
    //! since the number of nodes is even, zero is not an abscissa.
    //!
    //! declare
    //      external         f
    //      integer          j,n
    //      double precision a,b,result,drdeta,drdtheta, &
    //                       drdeta2,drdtheta2,drdetadtheta,par(n), &
    //                       absc,wg(20),xg(20),fval,deta,dtheta, &
    //                       deta2,dtheta2,detadtheta
    //
    //
    //! the abscissae and weights are given for the interval (0,+inf).
    //! xg     - abscissae of the 20-point gauss-laguerre rule
    //! wg     - weights of the 20-point gauss rule. since f yet
    //!          includes the weight function, the values in wg
    //!          are actually exp(xg) times the standard
    //!          gauss-laguerre weights
    //!
    //! abscissae and weights were evaluated with 100 decimal digit arithmetic.
    //
    //      data xg (  1) /   7.05398896919887533666890045842150958d-2 /
    xg[1] = 0.07053988969198875;
    //      data xg (  2) /   3.72126818001611443794241388761146636d-1 /
    xg[2] = 0.37212681800161146;
    //      data xg (  3) /   9.16582102483273564667716277074183187d-1 /
    xg[3] = 0.9165821024832735;
    //      data xg (  4) /   1.70730653102834388068768966741305070d0 /
    xg[4] = 1.707306531028344;
    //      data xg (  5) /   2.74919925530943212964503046049481338d0 /
    xg[5] = 2.749199255309432;
    //      data xg (  6) /   4.04892531385088692237495336913333219d0 /
    xg[6] = 4.048925313850887;
    //      data xg (  7) /   5.61517497086161651410453988565189234d0 /
    xg[7] = 5.6151749708616165;
    //      data xg (  8) /   7.45901745367106330976886021837181759d0 /
    xg[8] = 7.459017453671064;
    //      data xg (  9) /   9.59439286958109677247367273428279837d0 /
    xg[9] = 9.594392869581096;
    //      data xg ( 10) /   1.20388025469643163096234092988655158d1 /
    xg[10] = 12.038802546964316;
    //      data xg ( 11) /   1.48142934426307399785126797100479756d1 /
    xg[11] = 14.81429344263074;
    //      data xg ( 12) /   1.79488955205193760173657909926125096d1 /
    xg[12] = 17.948895520519375;
    //      data xg ( 13) /   2.14787882402850109757351703695946692d1 /
    xg[13] = 21.47878824028501;
    //      data xg ( 14) /   2.54517027931869055035186774846415418d1 /
    xg[14] = 25.451702793186904;
    //      data xg ( 15) /   2.99325546317006120067136561351658232d1 /
    xg[15] = 29.93255463170061;
    //      data xg ( 16) /   3.50134342404790000062849359066881395d1 /
    xg[16] = 35.013434240479;
    //      data xg ( 17) /   4.08330570567285710620295677078075526d1 /
    xg[17] = 40.83305705672857;
    //      data xg ( 18) /   4.76199940473465021399416271528511211d1 /
    xg[18] = 47.6199940473465;
    //      data xg ( 19) /   5.58107957500638988907507734444972356d1 /
    xg[19] = 55.810795750063896;
    //      data xg ( 20) /   6.65244165256157538186403187914606659d1 /
    xg[20] = 66.52441652561575;
    //
    //      data wg (  1) /   1.81080062418989255451675405913110644d-1 /
    wg[1] = 0.18108006241898925;
    //      data wg (  2) /   4.22556767878563974520344172566458197d-1 /
    wg[2] = 0.422556767878564;
    //      data wg (  3) /   6.66909546701848150373482114992515927d-1 /
    wg[3] = 0.6669095467018481;
    //      data wg (  4) /   9.15352372783073672670604684771868067d-1 /
    wg[4] = 0.9153523727830737;
    //      data wg (  5) /   1.16953970719554597380147822239577476d0 /
    wg[5] = 1.169539707195546;
    //      data wg (  6) /   1.43135498592820598636844994891514331d0 /
    wg[6] = 1.431354985928206;
    //      data wg (  7) /   1.70298113798502272402533261633206720d0 /
    wg[7] = 1.7029811379850228;
    //      data wg (  8) /   1.98701589079274721410921839275129020d0 /
    wg[8] = 1.9870158907927473;
    //      data wg (  9) /   2.28663578125343078546222854681495651d0 /
    wg[9] = 2.286635781253431;
    //      data wg ( 10) /   2.60583472755383333269498950954033323d0 /
    wg[10] = 2.6058347275538334;
    //      data wg ( 11) /   2.94978373421395086600235416827285951d0 /
    wg[11] = 2.949783734213951;
    //      data wg ( 12) /   3.32539578200931955236951937421751118d0 /
    wg[12] = 3.3253957820093194;
    //      data wg ( 13) /   3.74225547058981092111707293265377811d0 /
    wg[13] = 3.742255470589811;
    //      data wg ( 14) /   4.21423671025188041986808063782478746d0 /
    wg[14] = 4.21423671025188;
    //      data wg ( 15) /   4.76251846149020929695292197839096371d0 /
    wg[15] = 4.76251846149021;
    //      data wg ( 16) /   5.42172604424557430380308297989981779d0 /
    wg[16] = 5.421726044245574;
    //      data wg ( 17) /   6.25401235693242129289518490300707542d0 /
    wg[17] = 6.254012356932421;
    //      data wg ( 18) /   7.38731438905443455194030019196464791d0 /
    wg[18] = 7.387314389054435;
    //      data wg ( 19) /   9.15132873098747960794348242552950528d0 /
    wg[19] = 9.15132873098748;
    //      data wg ( 20) /   1.28933886459399966710262871287485278d1 /
    wg[20] = 12.893388645939996;
    //
    //
    //!           list of major variables
    //!           -----------------------
    //!           absc   - abscissa
    //!           fval*  - function value
    //!           result - result of the 20-point gauss formula
    //
    //      result       = 0.0d0
    result = 0.0;
    //      drdeta       = 0.0d0
    drdeta = 0.0;
    //      drdtheta     = 0.0d0
    drdtheta = 0.0;
    //      drdeta2      = 0.0d0
    drdeta2 = 0.0;
    //      drdtheta2    = 0.0d0
    drdtheta2 = 0.0;
    //      drdetadtheta = 0.0d0
    drdetadtheta = 0.0;
    //
    //      do j=1,20
    for (j = 1; j <= 20; ++j) {
        //       absc = a+b*xg(j)
        absc = (a + (b * xg[j]));
        //
        //       call f(absc, par, n, fval, deta, dtheta, &
        //              deta2,dtheta2,detadtheta)
        f(absc, par, n, fval, deta, dtheta, deta2, dtheta2, detadtheta);
        //
        //       result       = result    + fval*wg(j)
        result = (result + (fval * wg[j]));
        //       drdeta       = drdeta    + deta*wg(j)
        drdeta = (drdeta + (deta * wg[j]));
        //       drdtheta     = drdtheta  + dtheta*wg(j)
        drdtheta = (drdtheta + (dtheta * wg[j]));
        //       drdeta2      = drdeta2   + deta2*wg(j)
        drdeta2 = (drdeta2 + (deta2 * wg[j]));
        //       drdtheta2    = drdtheta2 + dtheta2*wg(j)
        drdtheta2 = (drdtheta2 + (dtheta2 * wg[j]));
        //       drdetadtheta = drdetadtheta + detadtheta*wg(j)
        drdetadtheta = (drdetadtheta + (detadtheta * wg[j]));
        //      enddo
    }
    //
    //      result       = result*b
    result = (result * b);
    //      drdeta       = drdeta*b
    drdeta = (drdeta * b);
    //      drdtheta     = drdtheta*b
    drdtheta = (drdtheta * b);
    //      drdeta2      = drdeta2*b
    drdeta2 = (drdeta2 * b);
    //      drdtheta2    = drdtheta2*b
    drdtheta2 = (drdtheta2 * b);
    //      drdetadtheta = drdetadtheta*b
    drdetadtheta = (drdetadtheta * b);
    //      return
    return;
    //      end
    //
    //
}
} // namespace

DirectPoint direct(double den, double temp) {
    if (!(std::isfinite(den) && den > 0 && std::isfinite(temp) && temp > 0))
        throw std::invalid_argument("Invalid direct EOS density or temperature");
    DirectWork w;
    double ages{}, fk{}, dfk{};
    //        call etages(xni,zbar,temp,ages)
    etages(avo * den, 1.0, temp, ages);
    bool converged = false;
    //        do i=1,100
    for (int i = 0; i < 100; ++i) {
        //         call xneroot(0,den,temp,abar,zbar,ionized,potmult,ages,fk,dfk)
        w.xneroot(0, den, temp, 1.0, 1.0, 1, 0, ages, fk, dfk);
        //         if (dfk .eq. 0.0) goto 11
        if (!std::isfinite(dfk) || dfk == 0)
            break;
        //         ratio   = fk/dfk
        //         agesnew = ages - ratio
        //         z       = abs((agesnew - ages)/ages)
        //         ages    = agesnew
        const double ratio = fk / dfk, next = ages - ratio;
        const double change = std::abs(next - ages) / std::max(std::abs(ages), 1.0);
        ages = next;
        //         if (z .lt. eostol .or. abs(ratio) .le. fpmin) goto 20
        // Avoid the original division by zero at eta=0. The reference's
        // default eostol is 1e-13; use that tolerance, plus stagnation.
        if (change < 1e-13 || std::abs(ratio) <= 1e-14) {
            converged = true;
            break;
        }
    }
    if (!converged || !std::isfinite(ages))
        throw std::runtime_error("Timmes direct EOS chemical-potential solve failed");
    //         call xneroot(1,den,temp,abar,zbar,ionized,potmult,ages,fk,dfk)
    w.xneroot(1, den, temp, 1.0, 1.0, 1, 0, ages, fk, dfk);
    return {w.pep,
            w.eele + w.epos,
            w.sep,
            w.dpepdd,
            w.dpepdt,
            w.dpepddd,
            w.dpepddt,
            w.dpepdtt,
            w.dsepdt,
            w.dsepddd,
            w.etaele,
            w.detadd,
            w.detadt,
            w.detaddt,
            w.xnefer + w.xnpfer,
            w.dxneferdd + w.dxnpferdd,
            w.dxneferdt + w.dxnpferdt,
            w.dxneferddt + w.dxnpferddt};
}
} // namespace octotigerII::helmholtz::detail
