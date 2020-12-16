
import json 
import time
import rasterio

import sys
import math
import numpy
import rasterio
from rasterio import features

def distance(viewpoint, point_coordinate):
    dx = point_coordinate[0] - viewpoint[0]
    dy = point_coordinate[1] - viewpoint[1]
    dist = (dx**2 + dy**2)**0.5
    return dist

def slope(viewpoint, point, numpy_raster):
    point_coordinate = d.xy(point[0], point[1])
    dist = distance(viewpoint, point_coordinate)

    # get locations of points in grid
    vrow, vcol = d.index(viewpoint[0], viewpoint[1])
    
    # get heights at points 
    height_v = numpy_raster[vrow, vcol] + viewpoint[2]
    height_p = numpy_raster[point[0], point[1]]
    # get height differnce 
    dh = height_p - height_v
    return dh/dist

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

# REMOVE IN MY CODE
output_file = 'test_tas1.tif'




#-- read the needed parameters from the file 'params.json' (must be in same folder)
jparams = json.load(open('params.json'))
# jparams = json.load(open('params2.json'))

#-- load in memory the input grid
d = rasterio.open(jparams['input_file'])    

#-- fetch the viewpoints
viewpoints = []
for i,each in enumerate(jparams['viewpoints']):
    vp = (jparams['viewpoints'][i]['xy'][0], jparams['viewpoints'][i]['xy'][1], jparams['viewpoints'][i]['height'])
    viewpoints.append(vp)

#get distance from the view point that can be seen 
radius_view = jparams['maxdistance']

#get raster resolution
resolution = d.res[0]

# make list of row and column numbers 
col_nrs = list(range(d.shape[1]))
row_nrs = list(range(d.shape[0]))

# make nested list of raster points 
raster_points = [[x, y] for y in col_nrs for x in row_nrs]

#-- numpy of input
npi  = d.read(1)
#-- the results of the viewshed in npvs, all values=0
npvs = numpy.zeros(d.shape, dtype=numpy.int8)

# create dictionary to store lists 
view_dict = {}
i = 0
# get the viewpoint from the list one by one
for v in viewpoints:
    view_pixels = []
    #-- index of this point in the numpy raster
    vrow, vcol = d.index(v[0], v[1])
    #-- put that pixel with value 2
    npvs[vrow , vcol] = 2

    # loop over all the raster points 
    for point in raster_points:
        col = point[1]
        row = point[0]
        # get coordinate values for the point
        coordinate = d.xy(row, col)
        # calculate the distance between the point and the viewpoint
        dist = distance(v, coordinate)
        # check if the raster point is at a radius of the viewpoint  
        if dist <= (radius_view + 0.5 * resolution) and npvs[row , col] != 2 and dist >= (radius_view - 0.5 * resolution):
            npvs[row , col] = 5
        elif npvs[row , col] != 5 and npvs[row , col] != 2 and (dist > (radius_view + 0.5 * resolution) or dist < (radius_view - 0.5 * resolution)): 
            npvs[row , col] = 3

        # check for circle edge outside of the raster extend
        if (row == 0 or row == npi.shape[0]-1 or col == 0 or col == npi.shape[1]-1) and dist <= (radius_view + 0.5 * resolution):
            npvs[row , col] = 5
    
        # use bresenham function to find pixels in line from viewpoint to end
        if npvs[row , col] == 5: 
            view_pixels.append(Bresenham_with_rasterio(d, (vrow, vcol), point))
    
    # store list of pixels in a dictionary 
    view_dict[i] = view_pixels
    i += 1

for j in range(len(viewpoints)):
    nested_line_lst = view_dict[j]
    for line in nested_line_lst:
        tan_old = -99
        for pixel in line[1:]:
            tan_new = slope(viewpoints[j], pixel, npi)
            col = pixel[1]
            row = pixel[0]
            if tan_new >= tan_old and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 1
                tan_old = tan_new
            elif tan_new < tan_old and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 0


#-- write this to disk
with rasterio.open(output_file, 'w', 
                    driver='GTiff', 
                    height=npi.shape[0],
                    width=npi.shape[1], 
                    count=1, 
                    dtype=rasterio.uint8,
                    crs=d.crs, 
                    transform=d.transform) as dst:
    dst.write(npvs.astype(rasterio.uint8), 1)

print("Viewshed file written to '%s'" % output_file)