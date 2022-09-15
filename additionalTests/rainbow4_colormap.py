#!/usr/bin/env python3

# -----------------------------------------------
# A rainbow-like palatte similar to the one
# used in Iliev+2006 plots.
#
# usage: import this file into your script, and
# provide the kwarg cmap=rainbow4 to your plots.
#
# Originally written by Yves Revaz, adapted for
# this repository by Mladen Ivkovic
# -----------------------------------------------

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt


_r = [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12,
    26,
    40,
    55,
    69,
    82,
    95,
    109,
    123,
    138,
    152,
    165,
    179,
    184,
    191,
    197,
    203,
    208,
    214,
    220,
    225,
    230,
    236,
    242,
    248,
    255,
    254,
    254,
    254,
    254,
    253,
    253,
    253,
    253,
    252,
    252,
    252,
    252,
    252,
    252,
    252,
    252,
    252,
    253,
    253,
    253,
    253,
    254,
    254,
    254,
    254,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
]
_g = [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8,
    18,
    27,
    37,
    45,
    56,
    64,
    75,
    83,
    93,
    102,
    112,
    121,
    127,
    133,
    140,
    146,
    152,
    158,
    164,
    170,
    176,
    182,
    187,
    193,
    199,
    204,
    208,
    211,
    216,
    221,
    224,
    229,
    233,
    237,
    242,
    246,
    249,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    248,
    243,
    237,
    233,
    228,
    222,
    216,
    210,
    205,
    201,
    195,
    190,
    184,
    179,
    174,
    170,
    165,
    159,
    154,
    150,
    145,
    140,
    134,
    129,
    125,
    120,
    109,
    101,
    91,
    82,
    72,
    63,
    55,
    45,
    36,
    26,
    18,
    8,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8,
    18,
    26,
    36,
    45,
    55,
    63,
    72,
    82,
    91,
    101,
    109,
    120,
    123,
    128,
    133,
    138,
    142,
    146,
    152,
    155,
    160,
    165,
    170,
    174,
    179,
    190,
    201,
    211,
    222,
    233,
    243,
    255,
]
_b = [
    2,
    6,
    11,
    14,
    19,
    23,
    26,
    31,
    34,
    38,
    43,
    46,
    51,
    53,
    57,
    62,
    65,
    70,
    74,
    77,
    82,
    85,
    89,
    94,
    97,
    102,
    104,
    108,
    113,
    116,
    121,
    125,
    128,
    133,
    136,
    140,
    145,
    148,
    153,
    155,
    159,
    164,
    167,
    172,
    176,
    179,
    184,
    187,
    191,
    196,
    199,
    204,
    206,
    210,
    215,
    218,
    223,
    227,
    230,
    235,
    238,
    242,
    247,
    250,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    248,
    242,
    236,
    230,
    225,
    220,
    214,
    208,
    203,
    197,
    191,
    184,
    179,
    165,
    152,
    138,
    123,
    109,
    95,
    82,
    69,
    55,
    40,
    26,
    12,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12,
    26,
    40,
    55,
    69,
    82,
    95,
    109,
    123,
    138,
    152,
    165,
    179,
    184,
    191,
    197,
    203,
    208,
    214,
    220,
    225,
    230,
    236,
    242,
    248,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
    255,
]


def _get_palette():
    r = np.array(_r, dtype=float) / 255.0
    g = np.array(_g, dtype=float) / 255.0
    b = np.array(_b, dtype=float) / 255.0

    return r, g, b


def _get_colormap(name="rainbow4"):
    """
    return a matplolib color map from the palette
    """

    LUTSIZE = plt.rcParams["image.lut"]

    r, g, b = _get_palette()
    ncols = r.shape[0]

    red = []
    green = []
    blue = []
    alpha = []

    for i in range(ncols):
        x = i / float(ncols - 1)
        red.append((x, r[i], r[i]))
        green.append((x, g[i], g[i]))
        blue.append((x, b[i], b[i]))
        alpha.append((x, 0.85, 0.85))

    cmapdata = {"red": red, "green": green, "blue": blue, "alpha": alpha}
    cmap = LinearSegmentedColormap(name, cmapdata, LUTSIZE)

    return cmap


rainbow4 = _get_colormap("rainbow4")


if __name__ == "__main__":

    x = np.linspace(1e-6, 1.0, 512)
    y = np.linspace(1e-6, 1.0, 512)
    X, Y = np.meshgrid(x, y)
    Z = X ** 2 + Y ** 2

    plt.figure()
    plt.imshow(Z, cmap=rainbow4)
    plt.show()
