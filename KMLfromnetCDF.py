#!/usr/bin/env python3
# Contouring netCDF data to KML
# Created 9 August 2022 by Sam Gardner <stgardner4@tamu.edu>

from errno import EDEADLK
from ntpath import join
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
import copy

if __name__ == "__main__":
    ## PARAMETERS (edit these)
    inputFile = "hgt.mon.ltm.nc" # input file name
    levelsStep = 6 # spacing between contours. Must be greater than 1.
    colormap = "viridis" # matplotlib colormap (see https://matplotlib.org/stable/tutorials/colors/colormaps.html for options)
    title = "Janurary 850 hPa height climatology" # Title
    alpha = 0.5 # transparency of polygons in kml, between 0 and 1

    ## THE ACTUAL CODE (this shouldn't need to be edited under normal circumstances)
    outputFileName = inputFile.replace(".nc", "")
    dataset = xr.open_dataset(inputFile)
    dataset = dataset.isel(time=0)
    dataset = dataset.sel(level=850.0)
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
    filledAx.contour(lonData, dataset.lat, hgtData, levels=levels, zorder=2, transform=ccrs.PlateCarree(), colors="k", linewidths=0.5)
    polys = filledAx.contourf(lonData, dataset.lat, hgtData, levels=levels, zorder=1, transform=ccrs.PlateCarree(), cmap=colormap)
    filledAx.add_feature(cfeat.COASTLINE.with_scale("50m"))
    filledAx.set_axis_off()
    filledFig.set_size_inches(1920*px, 1080*px)
    filledFig.savefig("FILLED_"+outputFileName+".png")
    polyKml = simplekml.Kml()
    polyFig = plt.figure()
    polyAx = plt.axes(projection=ccrs.PlateCarree())
    polyAx.contour(lonData, dataset.lat, hgtData, levels=levels, zorder=2, transform=ccrs.PlateCarree(), colors="k", linewidths=0.5)
    j = 0
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
                kmlPoly = polyKml.newpolygon(name="Between "+str(level)+" m and "+str(levels[j+1]), outerboundaryis=edgesOfPoly[0])
                rgb = np.array(colormap(norm(level))[0:-1]) * 255
                rgb = [hex(int(np.round(value))).replace("0x", "") for value in rgb]
                rgb.append(hex(int(np.round(alpha * 255))).replace("0x", ""))
                rgb = "".join(reversed(rgb))
                kmlPoly.style.polystyle.color = rgb
                innerPath = edgesOfPoly[0]
                polyAx.plot(innerPath[:,0], innerPath[:,1], transform=ccrs.PlateCarree(), color="lime", linewidth=0.5)
            elif len(edgesOfPoly) >= 2:
                for k in range(len(edgesOfPoly)):
                    parentPath = mpath(edgesOfPoly[k])
                    childPaths = copy.deepcopy(edgesOfPoly)
                    childPaths.pop(k)
                    isValid = True
                    childPathsList = None
                    for childPath in childPaths:
                        if childPathsList == None:
                            childPathsList = [childPath]
                        else:
                            childPathsList = childPathsList + [childPath]
                        if not np.all(parentPath.contains_points(childPath)):
                            isValid = False
                        if isValid:
                            break
                    kmlPoly = polyKml.newpolygon(name="PARENT: Between "+str(level)+" m and "+str(levels[j+1]), outerboundaryis=list(map(tuple, parentPath.vertices)), innerboundaryis=childPathsList)
                    rgb = np.array(colormap(norm(level))[0:-1]) * 255
                    rgb = [hex(int(np.round(value))).replace("0x", "") for value in rgb]
                    rgb.append(hex(int(np.round(alpha * 255))).replace("0x", ""))
                    rgb = "".join(reversed(rgb))
                    kmlPoly.style.polystyle.color = rgb
                    polyAx.plot(parentPath.vertices[:,0], parentPath.vertices[:,1], transform=ccrs.PlateCarree(), color="red", linewidth=0.5)
                    for childPath in childPaths:
                        polyAx.plot(childPath[:,0], childPath[:,1], transform=ccrs.PlateCarree(), color="cyan", linewidth=0.5)
        j += 1
    polyAx.add_feature(cfeat.COASTLINE.with_scale("50m"))
    polyAx.set_extent([-180, 180, -90, 90])
    polyAx.set_axis_off()
    polyFig.set_size_inches(1920*px, 1080*px)
    polyFig.savefig("POLY_"+outputFileName+".png")
    with open("POLY_"+outputFileName+".kml", "w") as f:
        f.write(polyKml.kml())
