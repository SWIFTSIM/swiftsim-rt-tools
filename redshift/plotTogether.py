#!/usr/bin/env python3

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

plot_kwargs = {
    "facecolor": "red",
    "s": 16,
    "alpha": 0.6,
    "linewidth": 0.0,
    "marker": ".",
}

analytic_solution_kwargs = {"color": "grey", "lw": 4, "alpha": 0.6, "zorder":1}

def fit_func(x,a,b):
    return b * x**a


# File to read from
analyticfile = "analytic.dat"
finlatorfile = "finlator.dat"

# Functions to convert scale factors and redshift
a2z = lambda a: 1/a - 1
z2a = lambda z: 1/(z+1)


def extract_data(file):
    # Read in data
    data = np.loadtxt(file)
    start = 0

    step = data[start:,0]
    a = data[start:,1]
    photonEnergy = data[start:,2]
    volume = data[start:,3]
    BBtemp = data[start:,4]
    withShiftingBB = data[0,-1]

    # Calculate photon energy density
    photonEnergyDensity = photonEnergy / volume

    return a, photonEnergy, photonEnergyDensity


if __name__ in ("__main__"):
    a, finE, finErho = extract_data(finlatorfile)
    __, anE, anErho = extract_data(analyticfile)

    # Scaled expected solutions
    expected_E = fit_func(a, -1., 1.)
    expected_Erho = fit_func(a, -4., 1.)

    expected_E = expected_E / expected_E[0] * finE[0]
    expected_Erho = expected_Erho / expected_Erho[0] * finErho[0]

    # Create figure and axes
    fig = plt.figure(dpi=300, figsize=(10,5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # Plot photon energy
    ax1.plot(a, finE, label="Finlator (2009)", zorder=2) #, **plot_kwargs)
    ax1.plot(a, anE, label="Integration by parts", zorder=2, ls="--", c="k")
    ax1.plot(a, expected_E, label="Analytic $\propto a^{-1}$", **analytic_solution_kwargs)
    secax = ax1.secondary_xaxis("top", functions=(a2z, z2a))
    secax.set_xlabel("Redshift")

    ax1.yaxis.get_offset_text().set_position((-0.05, 1))
    ax1.set_xlabel("Scale factor")
    ax1.set_ylabel("Total photon energy [IU]")
    ax1.set_title("Total energy")
    ax1.legend()

    # Plot photon energy density
    ax2.plot(a, finErho, label="Finlator (2009)", zorder=2) #, **plot_kwargs)
    ax2.plot(a, anErho, label="Integration by parts", zorder=2, ls="--", c="k")
    ax2.plot(a, expected_Erho, label="Analytic $\propto a^{-4}$", **analytic_solution_kwargs)
    secax = ax2.secondary_xaxis("top", functions=(a2z, z2a))
    secax.set_xlabel("Redshift")

    ax2.yaxis.get_offset_text().set_position((-0.05, 1))
    ax2.set_xlabel("Scale factor")
    ax2.set_ylabel("Total photon energy density [IU]")
    ax2.set_title("Energy density")
    ax2.legend()


    # Set scientific notation on y-axes
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    
    # Clean up and save plot
    plt.tight_layout()
    plt.savefig("finlator_analytic_same_plot.png", dpi=300)

