import json 
import time
import rasterio

import sys
import math
import numpy
import rasterio
from rasterio import features

arr = numpy.array([[1, 2], [3, 4], [5,6]])
print(arr)
print('------------------------')
print(numpy.flip(arr,0))


'''
intersect = numpy.argwhere(re==1)
    if intersect[0] != viewpoint:
        np.flip(intersect,0)
'''

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
    #return re
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
    return outlist

def Bresenham_with_rasterio2(raster, viewpoint, end):
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


    

output_file = 'bresham_test.tif'
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

end = (420, 22)
print(end)
v1= viewpoints[0]
v = d.index(v1[0],v1[1])
print('viewpoint', v)


output_original = Bresenham_with_rasterio(d, v, end)
output = Bresenham_with_rasterio2(d, v, end) # my version


'''


x, y = ((numpy.where(output)))
XY = [i for i in zip(x,y)]
#print(XY)

print(type(XY))
print(XY[0])

if XY[0] in [(1,2),(340, 318)]:
    print('yes')


line = []

if v in XY:             # append viewpoint
    line.append(v)
    ind = XY.index(v)
    XY.pop(ind)

#t = (line[-1][0]+1 , line[-1][1])
#print(t)
print('line', line)

up = (line[-1][0]-1, line[-1][1])
right = (line[-1][0] , line[-1][1]+1)
down = (line[-1][0]+1 , line[-1][1])
left = (line[-1][0] , line[-1][1]-1)
ur = (line[-1][0]-1, line[-1][1]+1)     # up right
dr = (line[-1][0]+1, line[-1][1]+1)     # down right
dl = (line[-1][0]+1, line[-1][1]-1)     # down left
ul = (line[-1][0]-1, line[-1][1]-1)     # up left
print ('up', up)
print ('right', right)
print ('down', down)
print ('left', left)
print ('ur', ur)
print ('dr', dr)
print ('dl', dl)
print ('ul', ul)
#next_point = v




while len(XY) != 0:
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

#print(output_original)
#print("-----------------------------------------")
#print(line)
#print(line)
'''

if output_original == output:
    print('yey')





'''
intersect = numpy.argwhere(output==1)
#print(intersect)
#if intersect[0] != v:
 #   numpy.flip(intersect,0)

x, y = ((numpy.where(output)))
XY = [i for i in zip(x,y)]

print(XY)

#vrow, vcol = d.index(v1[0], v1[1])


with rasterio.open(output_file, 'w', 
                    driver='GTiff', 
                    height=npi.shape[0],
                    width=npi.shape[1], 
                    count=1, 
                    dtype=rasterio.uint8,
                    crs=d.crs, 
                    transform=d.transform) as dst:
    dst.write(intersect.astype(rasterio.uint8), 1)

'''