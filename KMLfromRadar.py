#!/usr/bin/env python3
# Create KML file from radar raw file
# Created 15 August 2022 by Sam Gardner <stgardner4@tamu.edu>


from urllib import request
import pyart
from os import path, listdir
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as mpcolors
from cartopy import crs as ccrs
import KMLfromnetCDF

if __name__ == "__main__":
    alpha = 1
    requestedField = "reflectivity"
    SQIMinimum = 0.38


    basePath = path.abspath(path.dirname(__file__))
    inputDir = path.join(basePath, "radarData")
    outputDir = path.join(basePath, "radarKML")
    for i in range(len(listdir(inputDir))):
        fileToRead = sorted(listdir(inputDir), reverse=True)[i]
        radar = pyart.io.read(path.join(inputDir, fileToRead)).extract_sweeps([0])
        if "normalized_coherent_power" in radar.fields.keys():
            if SQIMinimum > 0:
                sqiValid = radar.fields["normalized_coherent_power"]["data"]
                sqiValid = np.where(sqiValid > SQIMinimum, 0, 1)
                finalRefl = np.ma.masked_array(radar.fields[requestedField]["data"], mask=sqiValid)
                radar.add_field_like(requestedField, "req_filtered", finalRefl)
                despekFilter = pyart.correct.despeckle_field(radar, "req_filtered")
                requestedField = "req_filtered"
            else:
                despekFilter = None
        else:
            despekFilter = None
        requestedFieldData = radar.fields[requestedField]["data"]
        norm  = mpcolors.Normalize(vmin=5, vmax=75)
        cmap = plt.get_cmap("pyart_NWSRef")
        cmap.set_over("black")
        fig = plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree())
        rmd = pyart.graph.RadarMapDisplay(radar)
        rmd.plot_ppi_map(requestedField, cmap=cmap, norm=norm, title_flag=False, colorbar_flag=False, ax=ax, fig=fig, gatefilter=despekFilter, width=2*160*1000, height=2*160*1000, embellish=False)
        plotHandle = ax.get_children()[0]

        Re = 6378137 # Radius of the WGS84-projection globe, in meters
        xOfRdr = ((radar.longitude["data"][0]+360)*np.pi*Re*np.sin(np.deg2rad(90-radar.latitude["data"][0])))/180
        yOfRdr = (radar.latitude["data"][0]*np.pi*Re)/180
        radarRelX = plotHandle.get_coordinates()[:,:,0]
        radarRelY = plotHandle.get_coordinates()[:,:,1]
        absoluteX = radarRelX + xOfRdr
        absoluteY = radarRelY + yOfRdr
        gateLat = (absoluteY*180)/(np.pi*Re)
        gateLon = (absoluteX*180)/(np.pi*Re*np.cos(np.deg2rad(radar.latitude["data"][0])))
        gateLon[gateLon > 180] = gateLon[gateLon > 180] - 360


        expandedFieldData = np.zeros((radarRelX.shape))
        expandedFieldData[:-1, :-1] = requestedFieldData.filled(np.nan)
        expandedFieldData[-1, :] = expandedFieldData[0, :]
        expandedFieldData[:, -1] = expandedFieldData[:, 0]
        plt.close(fig)
        fig = plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree())
        contourfs = ax.contourf(gateLon, gateLat, expandedFieldData, cmap=cmap, norm=norm, levels=np.arange(5, 75, 5), transform=ccrs.PlateCarree())
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax.set_aspect(1)
        ax.set_position([0, 0, 1, 1])
        print(ax.get_extent())
        print((radar.longitude["data"][0]))
        print((radar.latitude["data"][0]))
        ax.set_extent([(radar.longitude["data"][0]-0.5), (radar.longitude["data"][0]+0.5), (radar.latitude["data"][0]-0.5), (radar.latitude["data"][0]+0.5)])
        fig.colorbar(contourfs)
        fig.savefig(f"RADIAL_{fileToRead}.png")
        KMLfromnetCDF.create_kml_from_contourfs(contourfs, cmap, norm, np.arange(5, 75, 5), 0.75, "dBZ", False, f"RADIAL_{fileToRead}")
