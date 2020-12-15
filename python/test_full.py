
import json 
import time
import rasterio

import sys
import math
import numpy
import rasterio
from rasterio import features

def Bresenham_with_rasterio(raster, start, end):
    d = raster
    a = start
    b = end 
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append(d.xy(a[0], a[1]))
    v["coordinates"].append(d.xy(b[0], b[1]))
    shapes = [(v, 1)]
    re = features.rasterize(shapes,
                            out_shape=d.shape,
                            all_touched=True,
                            transform=d.transform)
    out = numpy.argwhere(re==1)
    outlist = []
    for el in out:
        outlist.append(tuple(el))
    if a[0]>b[0] and a[1]<=b[1]:
        outlist = sorted(outlist, key=lambda x: (-x[0], x[1]))
    elif a[0]>b[0] and a[1]>=b[1]:
        outlist = sorted(outlist, key=lambda x: (-x[0], -x[1]))
    elif a[0]<=b[0] and a[1]>b[1]:
        outlist = sorted(outlist, key=lambda x: (x[0], -x[1]))
    print(outlist)
    return outlist

output_file = 'test_tas_both_fullcircles.tif'
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

# get the viewpoint from the list one by one
for v in viewpoints:

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
        # -----------   dist >= (radius_view - resolution) and 

        if dist <= (radius_view + resolution) and npvs[row , col] != 2:
            npvs[row , col] = 1
        elif npvs[row , col] != 1 and npvs[row , col] != 2: 
            npvs[row , col] = 3

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