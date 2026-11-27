import numpy as np
from astropy.io import fits
from scipy import ndimage
import os
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import subprocess
import re
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
import argparse

parser = argparse.ArgumentParser(
    prog="mos_region_select.py",
    description="Selects the source and background regions for mos instruments"
)

parser.add_argument("--interactive", action="store_true")
parser.add_argument("--not-interactive", dest="interactive", action="store_false")
parser.set_defaults(interactive=True)

parser.add_argument("-i", "--image", default="mos1-all.fits")
parser.add_argument("-e", "--exposure-map", default="mos1_expmap_2-10keV.fits")
parser.add_argument("--source-mode", choices=["m","a"], default="a")
parser.add_argument("--back-mode", choices=["m","a"], default="a")
parser.add_argument("--source-radius", type=int, default=30)
parser.add_argument("--back-radius", type=int, default=15)
parser.add_argument("--target-mode", choices=["m","a"], default="a")
parser.add_argument("--instrument", type=str, choices=["mos1", "mos2"], default="mos1")

args = parser.parse_args()

def source_selector(image, exposure, target_x, target_y,
    background_positions, background_radius,
    target_area, args=args):

    # intial  background totals, there will increase

    bkgd_cts = 0
    bkgd_exposure = 0

    # here we set arrays so every coord can be re;ated to anotyer

    yy, xx = np.indices(image.shape)

    # Calculate total counts and exposure in background circles

    for bx, by in background_positions:
        distance = np.hypot(xx - bx, yy - by)


        # this creates a masl that makrs which pixels should be included in bkgd calculation:

        # 1. only pixels in or on edge of blkgd circle are included
        #  2. pixels with valid values (no NaN or infinities)
        # 3. only pixels with positive exposure

        region = (
            (distance <= background_radius) &
            np.isfinite(exposure) &
            (exposure > 0)
        )

        # cumulative sum of all counts in background and exposure map

        bkgd_cts += np.sum(image[region])
        bkgd_exposure += np.sum(exposure[region])

    if bkgd_exposure <= 0:
        raise ValueError("No exposure in background regions")

    bkgd_rate = bkgd_cts / bkgd_exposure


    # determine the maximum radius that fits in the exposed area

    edge_distance = ndimage.distance_transform_edt(target_area)

    # round the coords 

    target_ix = int(round(target_x))
    target_iy = int(round(target_y))

    #  max radius of edge so it doesnt fall out of pixel

    max_radius = int(edge_distance[target_iy, target_ix])
    if max_radius < 1:
        raise ValueError("Target (centre of AGN) is too close to the edge")


    if args.interactive:
        mode = input(
            "Select source radius manually or automatically? (m/a): "
        ).strip().lower()
    else:
        mode = args.source_mode

    if mode == "m":
        if args.interactive:
            radius = float(input("Enter source radius in pixels: "))
        else:
            radius = args.radius

        if radius <= 0 or radius > max_radius:
            raise ValueError(
                f"Radius must be between 0 and {max_radius} pixels; the region must fall entirely withtin the area."
            )

    elif mode == "a":
        #where snr is the signal to noise ratio

        distances_to_background = [

            np.hypot(target_x - bx, target_y - by)

            for bx, by in background_positions

        ]

        max_radius_background = min(distances_to_background) - background_radius
        max_radius = min(max_radius, int(max_radius_background))

        if max_radius < 1:
            raise ValueError("No space for a source region without overlapping background.")

        radii = np.arange(1, max_radius + 1, 1)
        snr_values = []

        dist_from_target = np.hypot(
            xx - target_x,
            yy - target_y
        )

        for radius_test in radii:
            region = (
                (dist_from_target <= radius_test) &
                target_area &
                (exposure > 0)
            )
            source_cts = np.sum(image[region])
            source_exposure = np.sum(exposure[region])
            expected_background = bkgd_rate * source_exposure
            net_counts = source_cts - expected_background

            if source_cts + expected_background > 0:
                snr = net_counts / np.sqrt(
                    source_cts + expected_background
                )
            else:
                snr = 0

            snr_values.append(snr)

        best_index = np.argmax(snr_values)
        radius = radii[best_index]
        print("\nAutomatic source radius selected:")
        print("Radius:", radius, "pixels")
        print("Estimated signal/noise ratio:", round(snr_values[best_index], 2))

    else:
        raise ValueError("Please enter 'm' or 'a'.")

    return radius   

# Taking user input if interactive, taking from arguments if not
if args.interactive:
    imagefile = input ("Enter the name of your .fits image:")
else:
    imagefile = args.image

with fits.open(imagefile) as data: 
    image = data[0].data.astype(float)
    header = data[0].header

print("Image dimensions:", image.shape)
print("Max intensity:", image.max())
print("Min intensity:", image.min())

#these are coordinates that constrain our data analysis to the CCD 

#chip with the obs on it


#active_pixels = np.isfinite(image) & (image > 0)

if args.interactive:
    exposure_file = input("Enter the matching MOS1 exposure map filename: ")
else:
    exposure_file = args.exposure_map

with fits.open(exposure_file) as data:
    exposure = data[0].data.astype(float)
if exposure.shape != image.shape:
    raise ValueError("Image and exposure map are different shapes")

active_pixels = np.isfinite(exposure) & (exposure > 0)

print("There are,", np.sum(active_pixels), "active pixels")
print("There are,", np.sum(~active_pixels), "inactive pixels")

mean, median, std = sigma_clipped_stats(
    image,
    sigma=3.0,
    mask=~active_pixels
)

print("Estimated background:", median)
print("Estimated noise:", std)

image_subtracted = image - median

fwhm = 4
threshold = 5 * std

finder = DAOStarFinder(
    fwhm=fwhm,
    threshold=threshold,
)
source = finder(
    image_subtracted,
    mask=~active_pixels
)

if source is None or len(source) == 0:
    raise ValueError(
        "No sources detected - try lowering the detection threshold (enter nano mos1_region_select.py)"
    )
else:
    source.sort("flux")
    source = source[::-1][:20]
    print("y centroid:", source[0]["y_centroid"])
    print("Number of detected sources:", len(source))
    print("20 brightest sources:")
    print(source)


#develop a section to ask the user whether they want to select 

#the cetral target manualy or not

if args.interactive:
    selection_mode = input(
        "Select target automatically or manually? (a/m): "
    ).strip().lower()
else:
    selection_mode = args.target_mode

if selection_mode == "a":
    target = source[0]
    print("Automatically selected the brightest source.")


elif selection_mode == "m":
    target_id = int(input("Enter the ID of your target source: "))
    matches = source[source["id"] == target_id]

    if len(matches) == 0:
        raise ValueError("Target ID not found.")

    target = matches[0]
    print("Manually selected source.")

else:
    raise ValueError("Please enter 'a' or 'm'.")

target_x = float(target["x_centroid"])
target_y = float(target["y_centroid"])

print("Selected target coordinates:")
print("X:", target_x)
print("Y:", target_y)

#wcs is world coordinate system helps convert from image coords 

# to real units 

wcs = WCS(header)

#convert the targetx and targety (in pixels) to detector coords in x y

detx_degree, dety_degree = wcs.pixel_to_world_values(
    target_x, target_y
)

detx = detx_degree * 3600 / 0.05
dety = dety_degree * 3600 / 0.05

print("Detector x-coord:", detx)
print("Detector y-coord:", dety)

#now we use sas to convert the detector coords to sky coods 

command = [
    "ecoordconv",
    f"imageset={imagefile}",
    f"x={detx}",
    f"y={dety}",
    "coordtype=det",
    "withcoords=yes"
]

result = subprocess.run(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True
)

output = result.stdout + result.stderr

if result.returncode != 0:
    raise RuntimeError("ecoordconv failed:" + output)
output = result.stdout + result.stderr

if result.returncode != 0:
    raise RuntimeError("ecoordconv failed:" + output)

match = re.search(
    r"RA:\s*DEC:\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)",
    output
)
if match is None:
    raise RuntimeError("Could not find RA and DEC in SAS output:" + output)
ra = float(match.group(1))
dec = float(match.group(2))

print("Selected source sky coordinates:")
print("RA:", ra)
print("DEC:", dec)

#ensure that the bkgd and source are on same section of the decetor 

area_labels, n = ndimage.label(active_pixels)

target_ix = int(round(target_x))
target_iy = int(round(target_y))

target_label = area_labels[target_iy, target_ix]


if target_label == 0:
    raise ValueError("This source target is not on an exposed detector pixel, or may not be on same pixel as source")

print("Number of connected exposed areas:", n)
print("The target lies in the area:", target_label)

#now we will work on sleecting a bckgd based off of just pixels used in the "area" that the source lies within

if args.interactive:
    background_radius = int(
        input("Enter background circle radius in pixels: ")
    )
else:
    background_radius = args.back_radius

if background_radius <= 0:
    raise ValueError("Background radius must be greater than zero.")

target_exclusion = 35
source_exclusion = 15
min_separation = 25
#make sure that the right CCD is  being used, and that the bkgd region will be within that CCD onluy 

target_area = (area_labels == target_label)

valid_area = (
	target_area
	& np.isfinite(exposure)
	& (exposure>0)
)

yy, xx = np.ogrid[
    -background_radius:background_radius + 1,
    -background_radius:background_radius + 1
]

circle_structure = (
    xx**2 + yy**2 <= background_radius**2
)


possible_centres = ndimage.binary_erosion(
    valid_area,
    structure=circle_structure,
    border_value=0
)

candidate_y, candidate_x = np.where(possible_centres)

if len(candidate_x) == 0:
    raise ValueError(
        "No valid position for the background circle, select manually. "
    )

distances = np.hypot(
    candidate_x - target_x,
    candidate_y - target_y
)

best_index = np.argmax(distances)

background_x = int(candidate_x[best_index])
background_y = int(candidate_y[best_index])

selected_positions = [(background_x, background_y)]

print("\nSelected background circle:")
print("Centre:", background_x, background_y)
print("Radius:", background_radius, "pixels")
print("Distance from target:", round(distances[best_index], 2), "pixels")


source_radius = source_selector(
    image,
    exposure,
    target_x,
    target_y,
    selected_positions,
    background_radius,
    target_area
)

print("\nSelected source region:")
print("Centre:", target_x, target_y)
print("Radius:", source_radius, "pixels")

# just relaised i forgot to addd a dfunction to save rthe regions as .reg files

if args.interactive:
    instrument = input("Which instrument used? (mos1/mos2): ").strip().lower()
else:
    instrument = args.instrument

if instrument not in ("mos1", "mos2"):
    raise ValueError("Please enter mos1 or mos2.")
region_header = "# Region file format: DS9 version 4.1\nimage\n"

with open(f"{instrument}_source.reg", "w") as regfile:
    regfile.write(region_header)
    regfile.write(
        f"circle({target_x + 1:.3f},{target_y + 1:.3f},{source_radius:.3f})\n"
    )

with open(f"{instrument}_background.reg", "w") as regfile:
    regfile.write(region_header)
    regfile.write(
        f"circle({background_x + 1},{background_y + 1},{background_radius})\n"
    )

print(f"\nSaved {instrument}_source.reg")
print(f"Saved {instrument}_background.reg")



#Plotting the selected regions

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

fig, ax = plt.subplots(figsize=(10, 8))


ax.imshow(
    np.ma.masked_less_equal(image, 0),
    origin="lower",
    cmap="inferno",
    norm=LogNorm(vmin=1, vmax=float(image.max()))
)

# Mark the target source
ax.plot(
    target_x,
    target_y,
    marker="+",
    color="red",
    markersize=15,
    label="Target source"
)

source_circle = Circle(
    (target_x, target_y),
    source_radius,
    fill=False,
    edgecolor="red",
    linewidth=2,
    label="Source region"
)

ax.add_patch(source_circle)

# Draw the two background circles
for number, (x, y) in enumerate(selected_positions, start=1):

    circle = Circle(
        (x, y),
        background_radius,
        fill=False,
        edgecolor=["lime", "cyan"][number - 1],
        linewidth=2,
        label=f"Background {number}"
    )

    ax.add_patch(circle)
    ax.plot(x, y, "o", color=["lime", "cyan"][number - 1])


ax.set_xlabel("Image X (pixels)")
ax.set_ylabel("Image Y (pixels)")
ax.set_title("MOS1 image with selected background regions")
ax.legend()
# ax.set_xlim(530, 800)  
# ax.set_ylim(530, 800)
plt.tight_layout()

plt.savefig(f"{instrument}-regions.png")
# plt.show()
