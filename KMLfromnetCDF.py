#!/usr/bin/env python3
# Contouring netCDF data to KML
# Created 9 August 2022 by Sam Gardner <stgardner4@tamu.edu>

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path as mpath
from matplotlib import colors as pltcolors
import pandas as pd
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from cartopy.util import add_cyclic_point
import simplekml

def create_kml_from_contourfs(polys, cmap, norm, levels, alpha, units, showBorders, filename):
    debugFig = plt.figure()
    debugAx = debugFig.gca()
    polyKml = simplekml.Kml()
    j = 0
    parentPoly = None
    for polyCollection in polys.collections:
        level = levels[j]
        paths = polyCollection.get_paths()
        for path in paths:
            vertices = path.vertices
            codes = path.codes
            startIndices = np.where(codes == 1)[0]
            endIndices = np.where(codes == 79)[0]
            edgesOfPoly = list()
            for i in range(len(startIndices)):
                startIdx = startIndices[i]
                endIdx = endIndices[i]
                pathPoints = vertices[startIdx:endIdx+1]
                edgesOfPoly.append(pathPoints)
            if len(edgesOfPoly) == 1:
                kmlPoly = polyKml.newpolygon(name=f"Between {str(level)} {units} and {str(levels[j+1])} {units}", outerboundaryis=edgesOfPoly[0])
                rgb = np.array(cmap(norm(level))[0:-1]) * 255
                rgb = [hex(int(np.round(value))).replace("0x", "") for value in rgb]
                rgb.append(hex(int(np.round(alpha * 255))).replace("0x", ""))
                rgb = [rgbvalue.zfill(2) for rgbvalue in rgb]
                rgb = "".join(reversed(rgb))
                kmlPoly.style.polystyle.color = rgb
                if showBorders:
                    kmlPoly.style.linestyle.width = 1
                else:
                    kmlPoly.style.linestyle.width = 0
            elif len(edgesOfPoly) >= 2:
                parentPath = edgesOfPoly[0]
                childPaths = edgesOfPoly[1:]
                parentPathTupleList = list(map(tuple, parentPath))
                childPathTupleList = [list(map(tuple, childPath)) for childPath in childPaths]
                parentPoly = polyKml.newpolygon(name=f"Between {str(level)} {units} and {str(levels[j+1])} {units}", outerboundaryis=parentPathTupleList, innerboundaryis=childPathTupleList)
                rgb = np.array(cmap(norm(level))[0:-1]) * 255
                rgb = [hex(int(np.round(value))).replace("0x", "") for value in rgb]
                rgb.append(hex(int(np.round(alpha * 255))).replace("0x", ""))
                rgb = [rgbvalue.zfill(2) for rgbvalue in rgb]
                rgb = "".join(reversed(rgb))
                parentPoly.style.polystyle.color = rgb
                if showBorders:
                    kmlPoly.style.linestyle.width = 1
                else:
                    kmlPoly.style.linestyle.width = 0
                    kmlPoly.style.linestyle.color = "00000000"
        j += 1
    print("WRITING DATA...")
    with open(filename+".kml", "w") as f:
        f.write(polyKml.kml())

if __name__ == "__main__":
    ## PARAMETERS (edit these)
    inputFile = "hgt.mon.ltm.nc" # input file name
    levelsStep = 6 # spacing between contours. Must be greater than 1.
    colormap = "viridis" # matplotlib colormap (see https://matplotlib.org/stable/tutorials/colors/colormaps.html for options)
    title = "Janurary 850 hPa height climatology" # Title
    alpha = 0.5 # transparency of polygons in kml, between 0 and 1
    minLat = 1 # minimum latitude to plot
    maxLat = 90 # maximum latitude to plot
    minLon = 230 # minimum longitude to plot
    maxLon = 300 # maximum longitude to plot
    units = "m" # units of data

    ## THE ACTUAL CODE (this shouldn't need to be edited under normal circumstances)
    outputFileName = inputFile.replace(".nc", "")
    dataset = xr.open_dataset(inputFile)
    dataset = dataset.isel(time=0)
    dataset = dataset.sel(level=850.0)
    # dataset = dataset.sel(lat=slice(maxLat, minLat), lon=slice(minLon, maxLon))
    print("\n\n")
    print("DATASET VARIABLES: "+str(list(dataset.variables)))
    print("DATASET MIN: "+str(np.min(dataset.hgt.data)))
    print("DATASET MAX: "+str(np.max(dataset.hgt.data)))
    print("DATASET AVG: "+str(np.mean(dataset.hgt.data)))
    print(pd.DataFrame(dataset.hgt))
    startValue = np.round(np.min(dataset.hgt.data), decimals=0) + 1.0 - levelsStep
    finalValue = np.round(np.min(dataset.hgt.data), decimals=0) + 1.0 - levelsStep
    datasetMax = np.round(np.max(dataset.hgt.data), decimals=0) - 1.0 + levelsStep
    while finalValue <= datasetMax:
        finalValue = finalValue + levelsStep
    levels = np.arange(startValue, finalValue, levelsStep)
    colormap = plt.get_cmap(colormap)
    norm = pltcolors.Normalize(vmin=startValue, vmax=finalValue)
    linesFig = plt.figure()
    linesAx = plt.axes(projection=ccrs.PlateCarree())
    hgtData, lonData = add_cyclic_point(dataset.hgt.data, coord=dataset.lon.data)
    contours = linesAx.contour(lonData, dataset.lat, hgtData, levels=levels, zorder=1, transform=ccrs.PlateCarree(), colors="k", linewidths=0.5)
    linesAx.add_feature(cfeat.COASTLINE.with_scale("50m"))
    linesAx.set_axis_off()
    px = 1/plt.rcParams["figure.dpi"]
    linesFig.set_size_inches(1920*px, 1080*px)
    linesKml = simplekml.Kml()
    i = 0
    allPaths = dict()
    for contour in contours.collections:
        level = levels[i]
        paths = contour.get_paths()
        pathsAtLevel = list()
        for path in paths:
            pathsAtLevel.append(path)
            vertices = path.vertices
            allPaths[i] = pathsAtLevel
            linesKml.newlinestring(name=str(level)+" m", coords=[vertex for vertex in vertices])
            linesAx.scatter(vertices[:,0], vertices[:,1], s=1, color="lime", transform=ccrs.PlateCarree())
        i += 1
    linesFig.savefig("LINES_"+outputFileName+".png")
    with open("LINES_"+outputFileName+".kml", "w") as f:
        f.write(linesKml.kml())
    plt.close(linesFig)
    filledFig = plt.figure()
    filledAx = plt.axes(projection=ccrs.PlateCarree())
    cc = filledAx.contour(lonData, dataset.lat, hgtData, levels=levels, zorder=2, transform=ccrs.PlateCarree(), colors="k", linewidths=0.5)
    filledAx.clabel(cc, cc.levels, inline=True, fontsize=5)
    polys = filledAx.contourf(lonData, dataset.lat, hgtData, levels=levels, zorder=1, transform=ccrs.PlateCarree(), cmap=colormap)
    filledAx.add_feature(cfeat.COASTLINE.with_scale("50m"))
    filledAx.set_axis_off()
    filledFig.set_size_inches(1920*px, 1080*px)
    filledFig.savefig("FILLED_"+outputFileName+".png")
    create_kml_from_contourfs(polys, colormap, norm, levels, alpha, units, True, "POLY_"+outputFileName)
