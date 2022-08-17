/* (Physical) constants. */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>

/* Physical
 * ------------------------------- */

/* Hydrogen mass in g */
#define const_mh 1.67262171e-24
/* Electron mass in g */
#define const_me 9.10938e-28
/* erg/K */
#define const_kboltz 1.3806488e-16
/* erg * s */
#define const_planck_h 6.62606957e-27
/* cm / s */
#define const_speed_light_c 2.99792458e+10
/* year in s */
#define const_yr 31557600.
/* Solar luminosity in erg/s */
#define const_L_Sun 3.826e33
#define const_adiabatic_index 1.666667

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Code internals
 * ------------------------------- */
#define SUCCESS 1
#define TINY_NUMBER 1.0e-20

/* Extra check just in case some program already defined this macro */
#ifndef RT_NGROUPS
#define RT_NGROUPS 4
#endif

#endif /* define CONSTANTS_H */
