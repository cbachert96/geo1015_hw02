import json 
import time
import rasterio
from rasterio import features
import math
import sys
import numpy


def main():
    
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

    output_viewshed(d, viewpoints, jparams['maxdistance'], jparams['output_file'])


def output_viewshed(d, viewpoints, maxdistance, output_file):
    npi  = d.read(1)
    v = viewpoints[0]
    vrow, vcol = d.index(v[0], v[1])
    print(d.shape)

if __name__ == '__main__':
    main()