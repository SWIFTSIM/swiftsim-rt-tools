#!/usr/bin/env python3

# --------------------------------------
# Plot the ionization and recombination
# rates for various temperatures
# --------------------------------------

import numpy as np
from thermochemistry_rates import thermochemistry_rates
from matplotlib import pyplot as plt

T = np.logspace(0, 11, 10000)

# Recombination rate for H+ in units of cm^3 s^-1
A_Hp = thermochemistry_rates.A_Hp(T)
# Dielectronic recombination rate for He+ in units of cm^3 s^-1
A_d = thermochemistry_rates.A_d(T)
# Recombination rate for He+ in units of cm^3 s^-1
A_Hep = thermochemistry_rates.A_Hep(T)
# Recombination rate for He++ in units of cm^3 s^-1
A_Hepp = thermochemistry_rates.A_Hepp(T)
# collisional ionization rate for H0 in units of cm^3 s^-1
G_H0 = thermochemistry_rates.G_H0(T)
# collisional ionization rate for He0 in units of cm^3 s^-1
G_He0 = thermochemistry_rates.G_He0(T)
# collisional ionization rate for He+ in units of cm^3 s^-1
G_Hep = thermochemistry_rates.G_Hep(T)


fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax1.loglog(T, A_Hp, label="$\\alpha_{H^+}$")
ax1.loglog(T, A_d, label="$\\alpha_{d}$")
ax1.loglog(T, A_Hep, label="$\\alpha_{He^{+}}$")
ax1.loglog(T, A_Hepp, label="$\\alpha_{He^{++}}$")
ax1.loglog(T, G_H0, "--", label="$\\Gamma_{H^{0}}$")
ax1.loglog(T, G_He0, "--", label="$\\Gamma_{He^{0}}$")
ax1.loglog(T, G_Hep, "--", label="$\\Gamma_{He^{+}}$")
ax1.legend()
ax1.set_title("Rates ($\\alpha$: recombination. $\\Gamma$: collisional ioinization)")
ax1.set_xlabel("T [K]")
ax1.set_ylabel("rates [cm$^3$ s$^{-1}$]")
ax1.set_ylim(1e-14, 1e-7)

ax2 = fig.add_subplot(1, 2, 2)
ax2.loglog(
    T,
    A_Hp / (A_Hp + G_H0),
    label="$\\alpha_{H^{+}} / (\\alpha_{H^{+}} + \\Gamma_{H^{0}})$",
)
ax2.loglog(
    T[T > 1e3],
    ((A_Hep + A_d) / G_He0)[T > 1e3],
    label="$(\\alpha_{He^{+}} + \\alpha_{d}) / \\Gamma_{He^{0}}$ for T > 1000 K",
)
ax2max = (((A_Hep + A_d) / G_He0)[T > 1e3]).max()
ax2min = ((G_Hep / A_Hepp)[T > 1e3]).min()
ax2.loglog([5e3, 5e3], [ax2min, ax2max], "k", label="T = 5000K", zorder=-3)

ax2.loglog(T, G_Hep / A_Hepp, label="$(\\Gamma_{He^{+}} / \\alpha_{He^{++}})$")
ax2.legend()
ax2.set_xlabel("T [K]")
ax2.set_ylabel("coefficients")
ax2.set_ylim(ax2min, ax2max)

plt.show()
