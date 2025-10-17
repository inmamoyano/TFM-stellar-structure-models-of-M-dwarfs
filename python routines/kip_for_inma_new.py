"""
Plotting routines for Kippenhahn diagram plots
Need the routines kipp_data.py, mesa_data.py and mkipp.py in the same directory
"""

import matplotlib.pyplot as plt
import numpy as np
import mkipp

# plot of Temperature against time, independent decoration
# There are options for the code to include "decoration" in the plot (axis
# names, etc) but I don't like them so I do it on the side

# Initialize plot
fig = plt.figure()
axis = plt.gca()

# mkipp.kipp_plot returns an object containing
#   kipp_plot.contour_plot : the return value of matplotlibs contourf. Can be
#                            used to create a colorbar with plt.colorbar()
#   kipp_plot.histories    : list of history files read. data can be accesed from this
#                            using the get("column_name") function
#   kipp_plot.xlims        : limits of data in x coordinate
kipp_plot = mkipp.kipp_plot(
    mkipp.Kipp_Args(
        xaxis="star_age",           # can be star_age or model_number
        time_units="Myr",           # if star_age, it can be yr, Myr, Gyr
        logs_dirs=["results/mass0.6_alpha1.82_z0.0142/"], 
                  # path to your log directory that must contain history,
                  # profiles.index and profile files
        identifier="logT",   # any quantity that is inside the profile files 
                             # this will be the colour map
        log10_on_data=False, # set to True if you want to plot the log of the
                             # quantity above for the colour map
        levels=np.arange(0, 8, 0.05),  #the levels at which the colourmap will
                             # change colours (min, max, step for colour change)
        decorate_plot=False, # set to True if you want to use auto decoration
                             # (see details in the associated routines)
        save_file=False,     # I save the file separately
    ),
    axis=axis,
)

# Set up location and title for colourbar
cbar = plt.colorbar(kipp_plot.contour_plot, pad=0.05)
cbar.set_label("logT")

# Set up y axis
axis.set_ylabel("Mass (solar masses)")
#axis.set_yscale("log")
#axis.set_ylim((1e-8, 1))
axis.set_ylim((0, 0.6))


# Set up x axis
axis.set_xlabel("Star age (Myr)")
axis.set_xscale("log")
axis.set_xlim((1e-6, 1.4e4))

# Save figure into jpg with good resolution
plt.savefig("Kippenhahn_test.jpg", dpi=300)
