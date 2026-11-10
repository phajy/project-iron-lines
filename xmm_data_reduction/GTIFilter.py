from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Plot formatting
cm = 0.3937008
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["figure.figsize"] = (8*cm, 5*cm)
mpl.rcParams["lines.linewidth"] = 0.5
mpl.rcParams["font.size"] = 5
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["errorbar.capthick"] = 0.2
mpl.rcParams["errorbar.capsize"] = 1.5
mpl.rcParams["lines.markersize"] = 0.1

def FindBounds(data, sigma=1):

    mean = np.mean(data)
    std = np.std(data)

    expression = f"(RATE<={mean+std})&&(RATE>={mean-std})"

    return expression, mean, std

def PlotBounds(data, mean, std):

    fig, ax = plt.subplots()

    ax.set_xlabel("Time")
    ax.set_ylabel("Rate (counts s$^{-1}$)")

    ax.plot(data["TIME"], data["RATE"], c="k")
    ax.hlines(mean, min(data["TIME"]), max(data["TIME"]), color="r")
    ax.hlines([mean-std, mean+std], min(data["TIME"]), max(data["TIME"]), color="r", ls="--")

    cond = ((data["RATE"]<=mean+std) & (data["RATE"]>=mean-std))

    ax.plot(data["TIME"][cond], data["RATE"][cond], color="b")

    ax.set_xlim(min(data["TIME"]), max(data["TIME"]))

    return fig, ax

f = open("GTIFilters.txt", "w")

# ==========================================================================================

hdu = fits.open("pn_ltcrv.fits")

data = hdu[1].data

expression, mean, std = FindBounds(data["RATE"], 1)

PlotBounds(data, mean, std)

plt.savefig("pn_gti_filter.png")

f.write(expression+"\n")

# ==========================================================================================

hdu = fits.open("mos1_ltcrv.fits")

data = hdu[1].data

expression, mean, std = FindBounds(data["RATE"], 1)

PlotBounds(data, mean, std)

plt.savefig("mos1_gti_filter.png")

f.write(expression+"\n")

# ==========================================================================================

hdu = fits.open("mos2_ltcrv.fits")

data = hdu[1].data

expression, mean, std = FindBounds(data["RATE"], 1)

PlotBounds(data, mean, std)

plt.savefig("mos2_gti_filter.png")

f.write(expression+"\n")

f.close()