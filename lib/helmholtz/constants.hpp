// C++ translation of F. X. (Frank) Timmes's stellar EOS routines.
// Original scientific algorithms: Timmes & Arnett (1999), ApJS 125, 277;
// Timmes & Swesty (2000), ApJS 126, 501. See NOTICE.md for provenance.
// Original Fortran statements are retained immediately above their translations.
#pragma once
#include <cmath>
namespace octotigerII::helmholtz::detail {
//
//! mathematical and physical constants (in cgs,except e0 which is in ev)
//!
//! 2006 codata recommended values of the physical constants
//! by coehn & taylor
//
//
//! math constants
//      real*8, parameter:: pi       = 3.1415926535897932384d0, &
//                          eulercon = 0.577215664901532861d0, &
//                          a2rad    = pi/180.0d0,  rad2a = 180.0d0/pi
inline const double pi = 3.141592653589793;
inline const double eulercon = 0.5772156649015329;
inline const double a2rad = (pi / 180.0);
inline const double rad2a = (180.0 / pi);
//
//! physical constants
//      real*8, parameter :: g       = 6.6742867d-8, &
//                           h       = 6.6260689633d-27, &
//                           hbar    = 0.5d0 * h/pi, &
//                           qe      = 4.8032042712d-10, &
//                           avo     = 6.0221417930d23, &
//                           clight  = 2.99792458d10, &
//                           kerg    = 1.380650424d-16, &
//                           ev2erg  = 1.60217648740d-12, &
//                           kev     = kerg/ev2erg, &
//                           amu     = 1.66053878283d-24, &
//                           mn      = 1.67492721184d-24, &
//                           mp      = 1.67262163783d-24, &
//                           me      = 9.1093821545d-28, &
//                           rbohr   = hbar*hbar/(me * qe * qe), &
//                           fine    = qe*qe/(hbar*clight), &
//                           hion    = 13.605698140d0, &
//                           ssol    = 5.6704d-5, &
//                           asol    = 4.0d0 * ssol / clight, &
//                           weinlam = h*clight/(kerg * 4.965114232d0), &
//                           weinfre = 2.821439372d0*kerg/h, &
//                            rhonuc  = 2.342d14
inline const double g = 6.6742867e-08;
inline const double h = 6.6260689633e-27;
inline const double hbar = ((0.5 * h) / pi);
inline const double qe = 4.8032042712e-10;
inline const double avo = 6.022141793e+23;
inline const double clight = 29979245800.0;
inline const double kerg = 1.380650424e-16;
inline const double ev2erg = 1.6021764874e-12;
inline const double kev = (kerg / ev2erg);
inline const double amu = 1.66053878283e-24;
inline const double mn = 1.67492721184e-24;
inline const double mp = 1.67262163783e-24;
inline const double me = 9.1093821545e-28;
inline const double rbohr = ((hbar * hbar) / ((me * qe) * qe));
inline const double fine = ((qe * qe) / (hbar * clight));
inline const double hion = 13.60569814;
inline const double ssol = 5.6704e-05;
inline const double asol = ((4.0 * ssol) / clight);
inline const double weinlam = ((h * clight) / (kerg * 4.965114232));
inline const double weinfre = ((2.821439372 * kerg) / h);
inline const double rhonuc = 234200000000000.0;
//
//
//! astronomical constants
//      real*8, parameter :: msol    = 1.9892d33, &
//                           rsol    = 6.95997d10, &
//                           lsol    = 3.8268d33, &
//                           mearth  = 5.9764d27, &
//                           rearth  = 6.37d8, &
//                           ly      = 9.460528d17, &
//                           pc      = 3.261633d0 * ly, &
//                           au      = 1.495978921d13, &
//                           secyer  = 3.1558149984d7
inline const double msol = 1.9892e+33;
inline const double rsol = 69599700000.0;
inline const double lsol = 3.8268e+33;
inline const double mearth = 5.9764e+27;
inline const double rearth = 637000000.0;
inline const double ly = 9.460528e+17;
inline const double pc = (3.261633 * ly);
inline const double au = 14959789210000.0;
inline const double secyer = 31558149.984;
} // namespace octotigerII::helmholtz::detail
