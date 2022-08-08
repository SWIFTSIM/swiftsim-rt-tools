#!/usr/bin/env python3

# ----------------------------------------------
# Generate an output times list to your liking
# ----------------------------------------------

import yaml
import numpy as np

unit_l = None
unit_v = None
t_end = None
with open(r"ilievTest0part3.yml") as paramfile:
    params = yaml.load(paramfile, Loader=yaml.FullLoader)

    unit_l = params["InternalUnitSystem"]["UnitLength_in_cgs"]
    unit_l = float(unit_l)
    unit_v = params["InternalUnitSystem"]["UnitVelocity_in_cgs"]
    unit_v = float(unit_v)
    t_end = params["TimeIntegration"]["time_end"]
    t_end = float(t_end)
    dt_min = params["TimeIntegration"]["dt_min"]
    dt_min = float(dt_min)
    dt_max = params["TimeIntegration"]["dt_max"]
    dt_max = float(dt_max)


# Check that t_end is indeed 5.5 Myr
unit_t = unit_l / unit_v
unit_myr = unit_t / (3600 * 24 * 365 * 1e6)

if abs((t_end * unit_myr) / 5.5 - 1.0) > 0.001:
    print("tend is {0:.3e} Myr, should be 5.5!".format(t_end * unit_myr))
else:
    print("t_end is {0:.3g} Myr".format(t_end * unit_myr))


output_heating = (
    np.logspace(np.log10(dt_max * unit_myr * 2), np.log10(0.5), 100) / unit_myr
)
output_cooling = np.linspace(0.51, 5.5, 100) / unit_myr
# make sure we have no precision issues for final snap
output_cooling[-1] = t_end

# If you're rebuilding this output times list:
# I read this 'dt' out from a run, then used it to generate the output list
dt_heat = 8.577071e-02

current_t = 0.0
outputtimes = []
# first 20 snapshots
for i in range(20):
    current_t += dt_heat
    outputtimes.append(current_t)
for i in range(20):
    current_t += 2 * dt_heat
    outputtimes.append(current_t)
for i in range(20):
    current_t += 4 * dt_heat
    outputtimes.append(current_t)
for i in range(20):
    current_t += 8 * dt_heat
    outputtimes.append(current_t)
for i in range(20):
    current_t += 32 * dt_heat
    outputtimes.append(current_t)
# fill up until heating stops
while current_t * unit_myr < 0.5:
    current_t += 64 * dt_heat
    outputtimes.append(current_t)


# do 50 snapshots until 1 Myr, and 50 thereafter
dt_cool = 0.5 / unit_myr / 50
for i in range(50):
    current_t += dt_cool
    outputtimes.append(current_t)
dt_cool = 4.5 / unit_myr / 50
for i in range(50):
    current_t += dt_cool
    outputtimes.append(current_t)

if outputtimes[-1] > t_end:
    outputtimes[-1] = t_end

if outputtimes[-2] == t_end:
    print("error: second to last snapshot > time_end")
    quit()


with open(r"snaplist.txt", "w") as outfile:
    # THIS HEADER IS IMPORTANT
    outfile.write("# Time\n")
    for t in outputtimes:
        outfile.write("{0:.6g}\n".format(t))

print("written snaplist.txt")
