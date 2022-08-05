#!/usr/bin/env python3


import xarray as xr
from simplekml import Kml, OverlayXY, ScreenXY, Units, RotationXY, AltitudeMode, Camera
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, **kw):
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""

    kml = Kml()
    altitude = kw.pop('altitude', 2e7)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = kw.pop('name', 'overlay')
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)

def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    fig = plt.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax

if __name__ == "__main__":

    # PARAMETERS (edit these)
    inputFile = "hgt.mon.ltm.nc" # input file name
    levels = np.arange(1150, 1580, 6) # contour levels, formatted as min, max, step
    colormap = "rainbow" # matplotlib colormap (see https://matplotlib.org/stable/tutorials/colors/colormaps.html for options)
    title = "Janurary 850 hPa height climatology" # Title


    outputFileName = inputFile.replace("nc", "")
    dataset = xr.open_dataset(inputFile)
    dataset = dataset.isel(time=0)
    dataset = dataset.sel(level=850.0)
    print("\n\n")
    print("DATASET VARIABLES: "+str(list(dataset.variables)))
    print("DATASET MIN: "+str(np.min(dataset.hgt.data)))
    print("DATASET MAX: "+str(np.max(dataset.hgt.data)))
    print("DATASET AVG: "+str(np.mean(dataset.hgt.data)))
    print(pd.DataFrame(dataset.hgt))
    print("\n\n")

    fig, ax = gearth_fig(llcrnrlon=0, llcrnrlat=-90, urcrnrlon=359.9999999999, urcrnrlat=90, pixels=10240)
    ax.contourf(dataset.lon, dataset.lat, dataset.hgt, levels=levels)
    ax.set_axis_off()
    fig.savefig(outputFileName+"png")
    make_kml(llcrnrlon=0, llcrnrlat=-90, urcrnrlon=359.9999999999, urcrnrlat=90, figs=[outputFileName+"png"], colorbar=None, kmzfile=outputFileName+"kmz", name=title)