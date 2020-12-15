
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
    out = numpy.argwhere(re==1)
    outlist = []
    for el in out:
        outlist.append(tuple(el))
    if viewpoint[0]>end[0] and viewpoint[1]<=end[1]:
        outlist = sorted(outlist, key=lambda x: (-x[0], x[1]))
    elif viewpoint[0]>end[0] and viewpoint[1]>=end[1]:
        outlist = sorted(outlist, key=lambda x: (-x[0], -x[1]))
    elif viewpoint[0]<=end[0] and viewpoint[1]>end[1]:
        outlist = sorted(outlist, key=lambda x: (x[0], -x[1]))
    print(outlist[0])
    return outlist

output_file = 'test_tas_both_circles.tif'
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
row_nr = list(range(d.shape[0]))

# make nested list of raster points 
raster_points = [[x, y] for y in row_nr for x in col_nrs]

#-- numpy of input
npi  = d.read(1)
#-- the results of the viewshed in npvs, all values=0
npvs = numpy.zeros(d.shape, dtype=numpy.int8)

# create dictionary to store lists 
view_dict = {}
i = 1
# get the viewpoint from the list one by one
for v in viewpoints:
    view_pixels = []
    #-- index of this point in the numpy raster
    vrow, vcol = d.index(v[0], v[1])
    #-- put that pixel with value 2
    npvs[vrow , vcol] = 2

    # loop over all the raster points 
    for point in raster_points:
        col = point[0]
        row = point[1]
        # get coordinate values for the point
        coordinate = d.xy(row, col)
        # calculate the distance between the point and the viewpoint
        dx = coordinate[0] - v[0]
        dy = coordinate[1] - v[1]
        dist = (dx**2 + dy**2)**0.5
        # check if the raster point is at a radius of the viewpoint  
        if dist <= (radius_view + 0.5 * resolution) and npvs[row , col] != 2 and dist >= (radius_view - 0.5 * resolution):
            npvs[row , col] = 1
        elif npvs[row , col] != 1 and npvs[row , col] != 2 and (dist > (radius_view + 0.5 * resolution) or dist < (radius_view - 0.5 * resolution)): 
            npvs[row , col] = 3

        # use bresenham function to find pixels in line from viewpoint to end
        if npvs[row , col] == 1: 
            view_pixels.append(Bresenham_with_rasterio(d, v, point))
    # store list of pixels in a dictionary 
    view_dict[i] = view_pixels
    i += 1

print()
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