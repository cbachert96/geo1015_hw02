#-- my_code_hw02.py
#-- Assignment 02 GEO1015.2020
#-- Louise Spekking 
#-- 4256778 
#-- Carolin Bachert
#-- 5382998
#-- Maundri Prihanggo 
#-- 5151279 

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
    
    x_index, y_index = ((numpy.where(re)))
    index_l = [i for i in zip(x_index,y_index)]             # creates list of index tuples, but not in right order

    line = []

    line.append(viewpoint)                  # append viewpoint as start point
    ind = index_l.index(viewpoint)               # find index in list
    index_l.pop(ind)                             # remove item from list

    while len(index_l) != 0:                     # like a walk algorythm, finds the next adjecant cell, if there is none, look for diagonals
        up = (line[-1][0]-1, line[-1][1])
        right = (line[-1][0] , line[-1][1]+1)
        down = (line[-1][0]+1 , line[-1][1])
        left = (line[-1][0] , line[-1][1]-1)
        ur = (line[-1][0]-1, line[-1][1]+1)     # up right
        dr = (line[-1][0]+1, line[-1][1]+1)     # down right
        dl = (line[-1][0]+1, line[-1][1]-1)     # down left
        ul = (line[-1][0]-1, line[-1][1]-1)     # up left
        if up in index_l:
            line.append(up)
            ind = index_l.index(up)
            index_l.pop(ind)
        elif right in index_l:
            line.append(right)
            ind = index_l.index(right)
            index_l.pop(ind)
        elif down in index_l:
            line.append(down)
            ind = index_l.index(down)
            index_l.pop(ind)
        elif left in index_l:
            line.append(left)
            ind = index_l.index(left)
            index_l.pop(ind)
        elif ur in index_l:
            line.append(ur)
            ind = index_l.index(ur)
            index_l.pop(ind)
        elif dr in index_l:
            line.append(dr)
            ind = index_l.index(dr)
            index_l.pop(ind)
        elif dl in index_l:
            line.append(dl)
            ind = index_l.index(dl)
            index_l.pop(ind)
        elif ul in index_l:
            line.append(ul)
            ind = index_l.index(ul)
            index_l.pop(ind)

    return line


def output_viewshed(d, viewpoints, maxdistance, output_file):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster
     
    Input:
        d:            the input datasets (rasterio format)  
        viewpoints:   a list of the viewpoints (x, y, height)
        maxdistance:  max distance one can see
        output_file:  path of the file to write as output
        
    Output:
        none (but output GeoTIFF file written to 'output-file')
    """  
        #get distance from the view point that can be seen 
    #radius_view = jparams['maxdistance']

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
            if dist <= (maxdistance + 0.5 * resolution) and npvs[row , col] != 2 and dist >= (maxdistance - 0.5 * resolution):
                npvs[row , col] = 5
            elif npvs[row , col] != 5 and npvs[row , col] != 2 and (dist > (maxdistance + 0.5 * resolution) or dist < (maxdistance - 0.5 * resolution)): 
                npvs[row , col] = 3

            # check for circle edge outside of the raster extend
            if (row == 0 or row == npi.shape[0]-1 or col == 0 or col == npi.shape[1]-1) and dist <= (maxdistance + 0.5 * resolution):
                npvs[row , col] = 5
        
            # use bresenham function to find pixels in line from viewpoint to end
            if npvs[row , col] == 5: 
                view_pixels.append(Bresenham_with_rasterio(d, (vrow, vcol), point))
        
        # store list of pixels in a dictionary 
        view_dict[i] = view_pixels
        i += 1

    for j in range(len(viewpoints)):
        nested_line_lst = view_dict[j]
        v = viewpoints[j]

        # get height of viewpoint
        height_v_plus = v[2]
        # get location of viewpoint in raster
        vrow, vcol = d.index(v[0], v[1])
        # get height viewpoint location on raster
        height_v_raster = npi[vrow, vcol]
        # height if viewpoint total 
        height_v = height_v_plus + height_v_raster

        for line in nested_line_lst:
            # calcualte the height of the first point in the line next to the viewpoint 
            height_p = npi[line[1][0], line[1][1]]
            # get x y coordinates of the first point in the line next to the viewpoint 
            x, y = d.xy(line[1][0], line[1][1])

            # calculate initial height difference 
            dh = height_p - height_v
            # calculate initial distance between viewpoint and second point
            dist = distance(v, (x, y))
            # calcualte inital slope 
            slope = dh/dist

            # implying that the first pixel is always visible
            npvs[line[1][0], line[1][1]] = 1
            # get all pixels in the list exept for the viewpoint, the first point
            for pixel in line[2:]:
                x_pixel, y_pixel = d.xy(pixel[0], pixel[1])

                # calculate the new dx 
                dx = distance(v, (x_pixel, y_pixel))
                # use y = ax + b to get predicted hight
                current_h = slope * dx + height_v

                if npvs[pixel[0], pixel[1]] != 1 and npvs[pixel[0], pixel[1]] != 2:
                    if current_h < npi[pixel[0], pixel[1]]: # changed <= to < 
                        npvs[pixel[0], pixel[1]] = 1
                        dh = npi[pixel[0], pixel[1]] - height_v
                        dx = distance(v, (x_pixel, y_pixel))
                        slope = dh / dx
                    else: 
                        if npvs[pixel[0], pixel[1]] != 1 and npvs[pixel[0], pixel[1]] != 2:
                            npvs[pixel[0], pixel[1]] = 0

    # wirte to file 
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



# def Bresenham_with_rasterio():
#     # d = rasterio dataset as above
#     a = (10, 10)
#     b = (100, 50)
#     #-- create in-memory a simple GeoJSON LineString
#     v = {}
#     v["type"] = "LineString"
#     v["coordinates"] = []
#     v["coordinates"].append(d.xy(a[0], a[1]))
#     v["coordinates"].append(d.xy(b[0], b[1]))
#     shapes = [(v, 1)]
#     re = features.rasterize(shapes, 
#                             out_shape=d.shape, 
#                             # all_touched=True,
#                             transform=d.transform)
#     # re is a numpy with d.shape where the line is rasterised (values != 0)



