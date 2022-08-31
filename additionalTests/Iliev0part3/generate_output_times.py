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


# If you're rebuilding this output times list:
# I read this 'dt' out from a run, then used it to generate the output list
#  dt_heat = 1.597603e-10 # for min step
dt_heat = 3.994010e-11

current_t = 0.0
new_t = 0.0
outputtimes = []
step = 32

count = 0
while new_t < t_end:
    for i in range(12):
        #  print(current_t, step, dt_heat, step * dt_heat)
        new_t = current_t + step * dt_heat
        if new_t > t_end:
            break
        current_t = new_t
        outputtimes.append(current_t)
        count += 1
    step *= 2

# add final time
if outputtimes[-1] < t_end * (1.0 - 1e-3):
    outputtimes.append(t_end * (1.0 - 1e-3))


with open(r"snaplist.txt", "w") as outfile:
    # THIS HEADER IS IMPORTANT
    outfile.write("# Time\n")
    for t in outputtimes:
        outfile.write("{0:.6g}\n".format(t))

print("written snaplist.txt")
