#!/usr/bin/env python3

# ----------------------------------------------------
# Find the ionization equilibrium for a given
# internal energy of the gas.
#
# The issue is that the ionization state determines
# the mean molecular weight of the gas, while it
# depends on the gas temperature. Conversely, the
# temperature determines the ionization state of the
# gas if it is in ionization equilibrium. So we need
# to iterate!
# ----------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import unyt

# globally defined constants
import constants
from gas_functions import *

# Set some initial conditions
# --------------------------------

XH = 0.76  # mass fraction of all hydrogen (HI + HII)
XHe = 1.0 - XH  # mass fraction of all helium (HeI + HeII + HeIII)
epsilon = 1e-3  # convergence criterion: relative to expected internal energy
iter_max = 100  # max iterations for Newton-Raphson
verbose = False  # am I talkative?
npoints = 100  # how many internal energies between u_min and u_max to compute for?
u_min = 1e11 * unyt.erg / unyt.g  # minimal specific internal energy to use
u_max = 1e17 * unyt.erg / unyt.g  # maximal specific internal energy to use
# NOTE: We use specific internal energy, which
# has units of energy per unit mass. ALSO: don't
# forget to set the unyts!


plotkwargs = {"alpha": 0.5}


def newton_raphson_iteration(u_expect, u_guess, T_guess, mu_guess):
    """
    One iteration of the Newton-Raphson root 
    finding algorithm to find the correct values for
    temperature and ionization equilibrium by checking
    whether the obtained internal energy given by the
    temperature guess is close enough to the expected
    internal energy.

    u_expect: The correct (given) internal energy of the gas
    u_guess: The current guess for internal energy of the gas
    T_guess: current guess for temperature

    returs: T_next
        next temperature to iterate over
    """
    T_guess.convert_to_cgs()  # better safe than sorry

    if T_guess < 0 * unyt.K:
        T_guess = 0.1 * unyt.K
        print(" Warning: Got negative temperature, resetting.")

    # find next temperature guess by solving linear equation
    # m * T_next + n = u_expect - u_guess_new ~ 0
    # NOTE: we pretend that the function that we're looking the
    # root of is f(T) = u_expect - u_guess(T), with u_expect = const.
    # so df/dT = - du_guess/dT, therefore we add a minus sign here.
    m = -internal_energy_derivative(T_guess, mu_guess)
    n = u_expect - u_guess - m * T_guess
    T_next = -n / m

    return T_next


def find_temperature_given_internal_energy(u_expect):
    """
    Find the gas temperature for a given
    internal energy u_expect.
    """

    # get first estimate for temperature.
    # First assume we're fully neutral.
    XH0 = XH
    XHp = 0.0
    XHe0 = XHe
    XHep = 0.0
    XHepp = 0.0
    mu_guess = mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp)
    T_guess = gas_temperature(u_expect, mu_guess)

    # If we're above the temperature threshold with this guess,
    # assume we're fully ionized as first guess instead.
    if T_guess > constants.T_thresh:
        XH0 = 0.0
        XHp = XH
        XHe0 = 0.0
        XHep = 0.0
        XHepp = XHe
        mu_guess = mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp)
        T_guess = gas_temperature(u_expect, mu_guess)

    # get updated mean molecular weight
    XH0, XHp, XHe0, XHep, XHepp = get_mass_fractions(T_guess, XH, XHe)
    mu_guess = mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp)

    # get first internal energy guess
    u_guess = internal_energy(T_guess, mu_guess)

    niter = 0
    du = u_expect - u_guess
    du_old = u_expect - u_guess

    # start iteration
    while abs(du) >= epsilon * u_expect:
        niter += 1

        if niter > iter_max:
            print("Error: Iteration didn't converge")
            print("     u              = ", u_guess)
            print("     T              = ", T_guess)
            print("     u_expect - u   = ", du)
            print("     1 - u/u_expect = ", 1.0 - u_guess / u_expect)
            return T_guess

        if verbose:
            print("i = ", niter)
            print("     u_expect           = ", u_expect)
            print("     u                  = ", u_guess)
            print("     T                  = ", T_guess)
            print("     delta u            = ", u_expect - u_guess)
            print("     delta u / u_expect = ", 1.0 - u_guess / u_expect)

        # do a Newton-Raphson iteration
        T_next = newton_raphson_iteration(u_expect, u_guess, T_guess, mu_guess)

        # Given the new temperature guess, compute the
        # expected mean molecular weight
        XH0, XHp, XHe0, XHep, XHepp = get_mass_fractions(T_next, XH, XHe)
        mu_next = mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp)

        # now given the new temperature and mass fraction guess, update the
        # expected gas internal energy
        u_next = internal_energy(T_next, mu_next)

        # save the old internal energy, and get the new one
        du_old = du
        du = u_expect - u_next

        # if we're oscillating between positive and negative values,
        # try a bisection to help out
        if du_old.v * du.v < 0.0:

            T_next = 0.5 * (T_guess + T_next)
            XH0, XHp, XHe0, XHep, XHepp = get_mass_fractions(T_next, XH, XHe)
            mu_next = mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp)
            u_next = internal_energy(T_next, mu_next)

        # reset what "current values" are and iterate again
        T_guess = T_next
        mu_guess = mu_next
        u_guess = u_next
        du = u_expect - u_next

    if verbose:
        print("Finished after", niter, "iterations; diff", 1.0 - du / u_expect)

    return T_guess, mu_guess, XH0, XHp, XHe0, XHep, XHepp


def plot_solution(u, T, mu, XH0, XHp, XHe0, XHep, XHepp, XH, XHe):
    """
    Plot the solutions.
    All arguments are expected to be arrays, except XH and XHe,
    which are scalars: total mass fraction of H and He, respectively.
    """

    # compute theorietical mean molecular weights,
    # internal energies, and mass fractions
    # ---------------------------------------------

    n = T.shape[0]
    u_theory = np.empty((n), dtype=float)
    mu_theory = np.empty((n), dtype=float)
    XH0_theory = np.empty((n), dtype=float)
    XHp_theory = np.empty((n), dtype=float)
    XHe0_theory = np.empty((n), dtype=float)
    XHep_theory = np.empty((n), dtype=float)
    XHepp_theory = np.empty((n), dtype=float)

    Tmin = T.v.min()
    Tmax = T.v.max()
    T_theory = np.logspace(np.log10(Tmin), np.log10(Tmax), npoints) * unyt.K

    for i, t in enumerate(T_theory):
        XH0c, XHpc, XHe0c, XHepc, XHeppc = get_mass_fractions(t, XH, XHe)
        muc = mean_molecular_weight(XH0c, XHpc, XHe0c, XHepc, XHeppc)
        uc = internal_energy(t, muc)

        u_theory[i] = uc
        mu_theory[i] = muc
        XH0_theory[i] = XH0c
        XHp_theory[i] = XHpc
        XHe0_theory[i] = XHe0c
        XHep_theory[i] = XHepc
        XHepp_theory[i] = XHeppc

    fig = plt.figure(figsize=(12, 6), dpi=300)
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    ax1.loglog(T, u, label="obtained results")
    ax1.loglog(T_theory, u_theory, label="analytical", ls=":")
    ax1.set_xlabel("gas temperature [K]")
    ax1.set_ylabel("specific internal energy [erg/g]")
    ax1.legend()
    #  ax1.set_title("obtained results")

    ax2.semilogx(u, mu, label="obtained results")
    ax2.semilogx(u_theory, mu_theory, label="expected results", ls=":")
    #  ax2.set_xlabel("gas temperature [K]")
    ax2.set_xlabel("specific internal energy [erg/g]")
    ax2.set_ylabel("mean molecular weight [1]")
    ax2.legend()
    #  ax2.set_title("Expected results")

    ax3.semilogx(
        u, XH0 + XHp + XHe0 + XHep + XHepp, "k", label="total", ls="-", **plotkwargs
    )
    ax3.semilogx(u, XH0, label="$H^0$", ls=":", **plotkwargs)
    ax3.semilogx(u, XHp, label="$H^+$", ls="-.", **plotkwargs)
    ax3.semilogx(u, XHe0, label="$He^0$", ls=":", **plotkwargs)
    ax3.semilogx(u, XHep, label="$He^+$", ls="-.", **plotkwargs)
    ax3.semilogx(u, XHepp, label="$He^{++}$", ls="--", **plotkwargs)
    ax3.legend()
    ax3.set_xlabel("specific internal energy [erg/g]")
    #  ax3.set_xlabel("gas temperature [K]")
    ax3.set_ylabel("gas mass fractions [1]")
    ax3.set_title("obtained results")

    ax4.semilogx(
        u,
        XH0_theory + XHp_theory + XHe0_theory + XHep_theory + XHepp_theory,
        "k",
        label="total",
        ls="-",
        **plotkwargs,
    )
    ax4.semilogx(u_theory, XH0_theory, label="$H^0$", ls=":", **plotkwargs)
    ax4.semilogx(u_theory, XHp_theory, label="$H^+$", ls="-.", **plotkwargs)
    ax4.semilogx(u_theory, XHe0_theory, label="$He^0$", ls=":", **plotkwargs)
    ax4.semilogx(u_theory, XHep_theory, label="$He^+$", ls="-.", **plotkwargs)
    ax4.semilogx(u_theory, XHepp_theory, label="$He^{++}$", ls="--", **plotkwargs)
    ax4.legend()
    #  ax4.set_xlabel("gas temperature [K]")
    ax4.set_xlabel("specific internal energy [erg/g]")
    ax4.set_ylabel("gas mass fractions [1]")
    ax4.set_title("Expected results")

    plt.tight_layout()
    #  plt.show()
    plt.savefig("ionization_equilibrium.png")

    return


if __name__ == "__main__":

    T_res = np.empty((npoints), dtype=float) * unyt.K
    mu_res = np.empty((npoints), dtype=float)
    XH0_res = np.empty((npoints), dtype=float)
    XHp_res = np.empty((npoints), dtype=float)
    XHe0_res = np.empty((npoints), dtype=float)
    XHep_res = np.empty((npoints), dtype=float)
    XHepp_res = np.empty((npoints), dtype=float)

    # set dummy internal energy of gas
    #  ugas = np.linspace(u_min, u_max, npoints)
    ugas = np.logspace(np.log10(u_min.v), np.log10(u_max.v), npoints) * u_min.units

    for i, u in enumerate(ugas):
        T, mu, XH0, XHp, XHe0, XHep, XHepp = find_temperature_given_internal_energy(u)
        T_res[i] = T
        mu_res[i] = mu
        XH0_res[i] = XH0
        XHp_res[i] = XHp
        XHe0_res[i] = XHe0
        XHep_res[i] = XHep
        XHepp_res[i] = XHepp

    plot_solution(
        ugas, T_res, mu_res, XH0_res, XHp_res, XHe0_res, XHep_res, XHepp_res, XH, XHe
    )
