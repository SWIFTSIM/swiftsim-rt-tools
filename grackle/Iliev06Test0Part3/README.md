Iliev 06 Test 0 (part 3)
--------------------------------

Heat up a gas parcel for 0.5 Myr with a constant flux, then let it cool.

NOTE: additionally requires GSL.

Reproducing Iliev+2006 (ui.adsabs.harvard.edu/abs/2006MNRAS.369.1625I/) Test 0 part 3.
Details:

- hydrogen gas only.
- the gas density distribution is fixed.
- original test in paper was in 128^3 cells.
- source of ionizing radiation has black-body spectrum with T_bb = 10^5K
- initially the gas is fully neutral and at 100K.
- Gas number density nH = 1 cm^-3
- Apply a photoionizing flux of F = 10^12 photons/s/cm^2 with a blackbody
  spectrum for 0.5 Myr. Turn it off afterwards, and let gas cool for 5 Myr.

  Rather than the actual flux, we inject the energy (density), using E = |F|/c.
  Converted to energy per unit time, the flux is
  A) for 3 photoionizing groups:

    Bin   0:  3.288e+15 -  5.945e+15 [Hz]  Luminosity/cm^2 = 1.350e+01 [erg/s/cm^2]
    Bin   1:  5.945e+15 -  1.316e+16 [Hz]  Luminosity/cm^2 = 2.779e+01 [erg/s/cm^2]
    Bin   2:  1.316e+16 -  5.879e+17 [Hz]  Luminosity/cm^2 = 6.152e+00 [erg/s/cm^2]


