
import json 
import time
import rasterio

import sys
import math
import numpy
import rasterio
from rasterio import features

def Bresenham_with_rasterio(raster, viewpoint, end):
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append(raster.xy(viewpoint[0], viewpoint[1]))
    v["coordinates"].append(raster.xy(end[0], end[1]))
    shapes = [(v, 1)]
    re = features.rasterize(shapes,
                            out_shape=raster.shape,
                            all_touched=True,
                            transform=raster.transform)
    
    x, y = ((numpy.where(re)))
    XY = [i for i in zip(x,y)]             # creates list of index tuples, but not in right order

    line = []

    line.append(viewpoint)                  # append viewpoint as start point
    ind = XY.index(viewpoint)
    XY.pop(ind)

    while len(XY) != 0:                     # like a walk algorythm, finds the next adjecant cell, if there is none, look for diagonals
        up = (line[-1][0]-1, line[-1][1])
        right = (line[-1][0] , line[-1][1]+1)
        down = (line[-1][0]+1 , line[-1][1])
        left = (line[-1][0] , line[-1][1]-1)
        ur = (line[-1][0]-1, line[-1][1]+1)     # up right
        dr = (line[-1][0]+1, line[-1][1]+1)     # down right
        dl = (line[-1][0]+1, line[-1][1]-1)     # down left
        ul = (line[-1][0]-1, line[-1][1]-1)     # up left
        if up in XY:
            line.append(up)
            ind = XY.index(up)
            XY.pop(ind)
        elif right in XY:
            line.append(right)
            ind = XY.index(right)
            XY.pop(ind)
        elif down in XY:
            line.append(down)
            ind = XY.index(down)
            XY.pop(ind)
        elif left in XY:
            line.append(left)
            ind = XY.index(left)
            XY.pop(ind)
        elif ur in XY:
            line.append(ur)
            ind = XY.index(ur)
            XY.pop(ind)
        elif dr in XY:
            line.append(dr)
            ind = XY.index(dr)
            XY.pop(ind)
        elif dl in XY:
            line.append(dl)
            ind = XY.index(dl)
            XY.pop(ind)
        elif ul in XY:
            line.append(ul)
            ind = XY.index(ul)
            XY.pop(ind)
        else:
            print('something went wrong')
            break

    return line