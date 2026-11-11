import sys

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as mpl

from astropy.wcs import WCS
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

# Plot formatting
cm = 0.3937008
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["figure.figsize"] = (8*cm, 5*cm)
mpl.rcParams["lines.linewidth"] = 0.2
mpl.rcParams["font.size"] = 5
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["errorbar.capthick"] = 0.2
mpl.rcParams["errorbar.capsize"] = 1.5
mpl.rcParams["lines.markersize"] = 0.1

def FindSource(hdu):
    """
    Find the centre of the brightest source in an image.

    Returns:
        source: centre coordinates in image coords
        factor: conversion factor from image to detector coords 
    """

    # Taking the data from the hdu
    data = hdu.data

    # Detecting sources
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=1*std)
    sources = daofind(data - median)

    # Taking the brightest source detected
    source = sources[sources["mag"] == sources["mag"].min()]
    source = np.array([source["x_centroid"].data[0], source["y_centroid"].data[0]])

    return source

def FindCorners(data):
    """
    Find the corners of the window in the image.

    Returns:
        points : a 2D array of coordinates
        rows   : The rows that contain a slice of the window
        columns: The columns that contain a slice of the window
    """

    # Finding the row and column indices that contain the window
    rows = np.invert((data == np.zeros(633)).all(axis=1)).nonzero()[0]
    columns = np.invert((data == np.zeros(633)).all(axis=0)).nonzero()[0]

    # Finding the coordinates
    p1 = np.array([np.where(data[rows[-1]]!=0)[0][-1], rows[-1]])
    p2 = np.array([columns[-1], np.where(data[:,columns[-1]]!=0)[0][0]])
    p3 = np.array([np.where(data[rows[0]]!=0)[0][0], rows[0]])
    p4 = np.array([columns[0], np.where(data[:,columns[0]]!=0)[0][-1]])

    return np.array([p1, p2, p3, p4]), rows, columns

def FindLineBisector(source, points):
    """
    Find parameters for the line between points, the perpendicular line from this to a source, and the intercept between the two
    """

    # Vector between the points
    cornerLineGrad, cornerLineInt = FindConnectingLine(points[0], points[1])

    # Line from source to line between points
    sourceLineInt = source[1] + source[0]/cornerLineGrad

    # Intercept between line between points and perpendicular bisector through source
    interceptX = -(cornerLineInt-sourceLineInt)/(cornerLineGrad+1/cornerLineGrad)
    intercept = [interceptX, cornerLineGrad*interceptX+cornerLineInt]

    return intercept, cornerLineGrad, cornerLineInt, sourceLineInt

def FindConnectingLine(point1, point2):

    diff = point1 - point2
    grad = diff[1]/diff[0]
    intercept = point1[1] - grad*point1[0]

    return grad, intercept

def FindBackgroundRegions(points, source, distances, margin1, margin2):

    # Finding the two corners furthest from the source
    index = np.argpartition(distances, len(distances)-2)[-2:]

    # Find the closest point to the source on the furthest edge of the window
    bisectorIntercept, cornerLineGrad, cornerLineInt, sourceLineInt = FindLineBisector(source, points[index])

    regions = []

    for idx in index:

        # Taking the selected and opposite corners
        point = points[idx]
        oppidx = np.array([idx - 2, idx + 2], dtype=int)[(np.array([idx - 2, idx + 2], dtype=int)<=3) & (np.array([idx - 2, idx + 2], dtype=int)>=0)][0]
        oppPoint = points[oppidx]

        # Finding the line connecting the two corners
        grad, intercept = FindConnectingLine(point, oppPoint)

        # The diameter of the region in the selected corner
        diameter = np.linalg.norm(point-bisectorIntercept) - margin1

        # Distance to travel along connecting line
        distance = np.sqrt(2*(diameter)**2)/2

        # Coefficients to solve for the x coordinate of the region centre
        a = 1+grad**2
        b = 2*grad*intercept - 2*grad*point[1] - 2*point[0]
        c = point[0]**2 + point[1]**2 + intercept**2 - 2*intercept*point[1] - distance**2

        # Possible solutions
        xs = np.array(((-b - np.sqrt(b**2 - 4*a*c)) / (2*a), (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)))

        # Only taking the solution that lies within the window
        x = xs[((xs<point[0]) & (xs>oppPoint[0])) | ((xs>point[0]) & (xs<oppPoint[0]))][0]
        
        # Taking the margin from the edge of the window
        diameter -= margin2

        # Storing the found region
        regions.append((x, grad*x+intercept, diameter/2))

    return regions

def FindRegions(hdu, margin1, margin2):
    """
    Find source and background regions within a small window EPIC-pn image

    Returns:
        srcRegion
        bgdRegions
        rows
        columns
    """

    # Finding the source
    source = FindSource(hdu)

    # Finding the corners of the image
    points, rows, columns = FindCorners(hdu.data)

    # Distances between every corner and the source
    distances = np.linalg.norm(np.full((4,2), source)-points, axis=1)

    # Finding the background regions
    bgdRegions = FindBackgroundRegions(points, source, distances, margin1, margin2)

    # Finding the distance from the source to the closest window edge
    sourceRegionIntercept, *_ = FindLineBisector(source, points[np.argpartition(distances, len(distances)-2)[:2]])
    radius = np.linalg.norm(sourceRegionIntercept-source) - margin2

    return [*source, radius], bgdRegions, rows, columns

def PlotRegions(hdu, srcRegion, bgdRegions, columns, rows):
    
    fig, ax = plt.subplots()

    # Plotting the image
    ax.imshow(hdu[0].data, origin='lower', norm=LogNorm(), cmap="plasma")

    # Cropping to the small window
    ax.set_xlim(columns[0]-2, columns[-1]+2)
    ax.set_ylim(rows[0]-2, rows[-1]+2)

    ax.set_xticks([])
    ax.set_yticks([])

    ax.add_patch(plt.Circle(srcRegion[0:2], srcRegion[2], color='lime', fill=False, lw=1))

    for region in bgdRegions:
        ax.add_patch(plt.Circle(region[0:2], region[2], color='lime', fill=False, lw=1))

    return fig, ax

# Opening the file
hdu = fits.open(f"{sys.argv[1]}-all.fits")

factor = hdu[0].header["CDELT1L"]

margin1 = float(sys.argv[2])
margin2 = float(sys.argv[3])
sourceMax = float(sys.argv[4])

srcRegion, bgdRegions, rows, columns = FindRegions(hdu[0], margin1, margin2)

if srcRegion[2]*factor > sourceMax:
    srcRegion[2] = sourceMax/factor

f = open(f"{sys.argv[1]}-regions.txt", "w")

f.write(f"{factor*srcRegion[0]},{factor*srcRegion[1]},{factor*srcRegion[2]}\n")
for region in bgdRegions:
    f.write(f"{factor*region[0]},{factor*region[1]},{factor*region[2]}\n")

# ==========================================================================================

PlotRegions(hdu, srcRegion, bgdRegions, columns, rows)

plt.savefig(f"{sys.argv[1]}-regions.png")