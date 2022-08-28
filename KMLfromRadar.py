#!/usr/bin/env python3
# Create KML file from radar raw file
# Created 15 August 2022 by Sam Gardner <stgardner4@tamu.edu>


import pyart
from os import path, listdir
from metpy.plots import ctables
import numpy as np
import simplekml
from datetime import datetime as dt, timedelta

if __name__ == "__main__":
    alpha = 1
    requestedField = "reflectivity"


    basePath = path.abspath(path.dirname(__file__))
    inputDir = path.join(basePath, "radarData")
    outputDir = path.join(basePath, "radarKML")
    for i in range(len(listdir(inputDir))):
        fileToRead = sorted(listdir(inputDir), reverse=True)[i]
        radar = pyart.io.read(path.join(inputDir, fileToRead)).extract_sweeps([0])
        requestedFieldData = radar.fields[requestedField]["data"]
        norm, cmap = ctables.registry.get_with_steps("NWSReflectivity", 5, 5)
        cmap.set_over("black")
        cmap.set_under("white")


        Re = 6378137 # Radius of the WGS84-projection globe, in meters
        ranges = radar.range["data"].reshape(1, -1)
        sampleLength = ranges[0,1]
        groundLength = sampleLength*np.cos(np.deg2rad(radar.elevation["data"].reshape(-1, 1)))
        rayWidthDeg = 360/radar.rays_per_sweep["data"][0]
        nearWidth = 2*(ranges-(groundLength/2))*np.tan(np.deg2rad(rayWidthDeg/2))
        farWidth = 2*(ranges+(groundLength/2))*np.tan(np.deg2rad(rayWidthDeg/2))
        xdisplacement = ranges*np.sin(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))
        ydisplacement = ranges*np.cos(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))
        xOfRdr = ((radar.longitude["data"][0]+360)*np.pi*Re*np.sin(np.deg2rad(90-radar.latitude["data"][0])))/180
        yOfRdr = (radar.latitude["data"][0]*np.pi*Re)/180
        
        nearWidthX = nearWidth*np.cos(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))
        farWidthX = farWidth*np.cos(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))
        groundLengthX = groundLength*np.sin(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))
        nearWidthY = nearWidth*np.sin(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))
        farWidthY = farWidth*np.sin(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))
        groundLengthY = groundLength*np.cos(np.deg2rad(radar.azimuth["data"].reshape(-1, 1)))

        radarRelX = np.zeros((requestedFieldData.transpose().shape[0] + 1, requestedFieldData.transpose().shape[1] + 1))
        radarRelY = np.zeros((requestedFieldData.transpose().shape[0] + 1, requestedFieldData.transpose().shape[1] + 1))
        testVar = True
        for rayNum in range(len(radar.azimuth["data"])):
            for gateNum in range(len(radar.range["data"])):
                radarRelX[gateNum, rayNum] = xdisplacement[rayNum, gateNum] - (nearWidthX[rayNum, gateNum]/2)
                radarRelY[gateNum, rayNum] = ydisplacement[rayNum, gateNum] - (groundLengthY[rayNum]/2)
            radarRelX[-1, rayNum] = xdisplacement[rayNum, gateNum] - (farWidthX[rayNum, gateNum]/2)
            radarRelY[-1, rayNum] = ydisplacement[rayNum, gateNum] + (groundLengthY[rayNum]/2)
        for gateNum in range(len(radar.range["data"])):
            radarRelX[gateNum, -1] = xdisplacement[-1, gateNum] + (nearWidthX[-1, gateNum]/2)
            radarRelY[gateNum, -1] = ydisplacement[-1, gateNum] - (groundLengthY[-1]/2)
        radarRelX[-1, -1] = xdisplacement[-1, -1] + (farWidthX[-1, -1]/2)
        radarRelY[-1, -1] = ydisplacement[-1, -1] + (groundLengthY[-1]/2)        
        
        absoluteX = radarRelX + xOfRdr
        absoluteY = radarRelY + yOfRdr
        gateLat = (absoluteY*180)/(np.pi*Re)
        gateLon = (absoluteX*180)/(np.pi*Re*np.cos(np.deg2rad(radar.latitude["data"][0])))
        
        radialKml = simplekml.Kml()
        startTime = dt.utcnow()
        for rayNum in range(requestedFieldData.shape[0]):
            print(f"Ray {rayNum} of {requestedFieldData.shape[0]}")
            print(f"Time: {dt.utcnow() - startTime}")
            startTime = dt.utcnow()
            for gateNum in range(requestedFieldData.shape[1]):
                level = requestedFieldData[rayNum, gateNum]
                if level is np.ma.masked:
                    continue
                edgesOfPoly = [(gateLat[gateNum, rayNum], gateLon[gateNum, rayNum]), # lower left
                                (gateLat[gateNum+1, rayNum], gateLon[gateNum+1, rayNum]), # upper left
                                (gateLat[gateNum+1, rayNum+1], gateLon[gateNum+1, rayNum+1]), # upper right
                                (gateLat[gateNum, rayNum+1], gateLon[gateNum, rayNum+1])] # lower right
                kmlPoly = radialKml.newpolygon(name=f"{str(level)} dBZ", outerboundaryis=edgesOfPoly)
                rgb = np.array(cmap(norm(level))[0:-1]) * 255
                rgb = [hex(int(np.round(value))).replace("0x", "") for value in rgb]
                rgb.append(hex(int(np.round(alpha * 255))).replace("0x", ""))
                rgb = "".join(reversed(rgb))
                kmlPoly.style.polystyle.color = rgb
        with open("RADIAL_"+fileToRead+".kml", "w") as f:
            f.write(radialKml.kml())