''''
        # find value of the absolute altitude of the viewpoint
        height_vp = npi[vrow, vcol] + v[2]
    
        # use bresenham function to find pixels in line from viewpoint to end
        if npvs[row , col] == 5: 
            line = Bresenham_with_rasterio(d, (vrow, vcol), point)
            
            old_slope = -99

            for pixel in line[1:]:
                new_slope = slope(v, pixel, npi)

                if new_slope > old_slope and npvs[row , col] != 2:
                    npvs[row , col] = 1
                    old_slope = new_slope
                
                elif new_slope < old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                    npvs[row , col] = 0
                
                elif new_slope == old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                    npvs[row , col] = 0

                    

    # store list of pixels in a dictionary 
    view_dict[i] = view_pixels
    i += 1

for j in range(len(viewpoints)):
    nested_line_lst = view_dict[j]
    for line in nested_line_lst:
        tan_old = -99
        for k in range(len(line)-1):
            tan_new = slope(viewpoints[j], line[k+1], npi)
            col = line[k+1][0]
            row = line[k+1][1]
            if tan_new >= tan_old and npvs[row , col] != 2:
                npvs[row , col] = 1
                tan_old = tan_new
            elif tan_new < tan_old and npvs[row , col] != 2:
                npvs[row , col] = 0
''' 


           old_slope = -99

            for pixel in line[1:]:
                new_slope = slope(v, pixel, npi)

                if new_slope > old_slope and npvs[row , col] != 2:
                    npvs[row , col] = 1
                    old_slope = new_slope
                
                elif new_slope < old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                    npvs[row , col] = 0
                
                elif new_slope == old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                    npvs[row , col] = 0




        old_slope = -99
        for pixel in line[1:]:
            new_slope = slope(v, pixel, npi)

            print('new slope', new_slope)

            if new_slope > old_slope and npvs[row , col] != 2:
                print('yes')
                npvs[row , col] = 1
                old_slope = new_slope
                print('old slope', old_slope)
            elif new_slope < old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 0
                print('no 1')

            elif new_slope == old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 0    
                print('no 2')



# store list of pixels in a dictionary 
    view_dict[i] = view_pixels
    i += 1

for j in range(len(viewpoints)):
    v_

    nested_line_lst = view_dict[j]
    for line in nested_line_lst:
        old_slope = -99
        for pixel in line[1:]:
            new_slope = slope(v, pixel, npi)

            print('new slope', new_slope)

            if new_slope > old_slope and npvs[row , col] != 2:
                print('yes')
                npvs[row , col] = 1
                old_slope = new_slope
                print('old slope', old_slope)
            elif new_slope < old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 0
                print('no 1')

            elif new_slope == old_slope and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 0    
                print('no 2')




for i in range(len(viewpoints)):
    v = viewpoints[i]
    #-- index of this point in the numpy raster
    vrow, vcol = d.index(v[0], v[1])

    # get height of viewpoint 
    v_height = npi[vrow, vcol] + v[2]

    # get nested list from dictionary
    view_pixel_lst = view_dict[i]

    for line in view_pixel_lst:
        dis_init = distance(v, line[1])
        slope_init = ((npi[line[1]])-v_height) / dis_init
        tcur = (slope_init * dis_init) + v_height

        for pixel in line[2:]:
            slope_new = (npi[pixel]- v_height) / distance(v,pixel)

            if slope_new > slope_init and npvs[row , col] != 2:
                print('yes')
                npvs[row , col] = 1
                slope_init = slope_new
                print('init slope', slope_init)
            elif slope_new < slope_init and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 0
                print('no 1')
            elif slope_new == slope_init and npvs[row , col] != 2 and npvs[row , col] != 1:
                npvs[row , col] = 0    
                print('no 2')