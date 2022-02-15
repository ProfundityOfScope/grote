#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 23:27:59 2020
a
@author: bruzewskis
"""

from astropy.table import Table, join
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

import schedmaster as sm
from optscheduler import write_OPT

if __name__=='__main__':
    
    # Prep data table
    colnames = ['name','priority','system','epoch','ra','dec', 
                'lsr','regime','vel','blank']
    datdict = { c:[] for c in colnames}
    with open('inputs/20A-439.pst', 'r') as fp:
        for line in fp:
            sline = line.strip().split(';')[:-1]
            for i in range(len(sline)):
                datdict[ colnames[i] ].append( sline[i] )
    data = Table(datdict)
    
    # Skim down a bit
    data = data[['name','ra','dec','priority']]
    data['name'] = [ n.replace('4FGL ','') for n in data['name'] ]
    data['ra'] = Angle(data['ra'], unit='h').deg * u.deg
    data['dec'] = Angle(data['dec'], unit='deg').deg * u.deg
    data['priority'] = [ p[-2] for p in data['priority'] ]
    
    # Read in calibrators
    cal_file = 'inputs/calcat_L.fits'
    cals = Table.read(cal_file)
    
    # set up skycoords for ease of use
    scs = SkyCoord(data['ra'], data['dec'])
    
    blocks = []
    for i in range(len(data)):
        block = sm.schedule_block(scs[i], data['name'][i], cals)
        blocks.append(block)
    write_OPT(blocks, scs, data['priority'])