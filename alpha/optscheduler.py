#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SimEVLA.py: This module simulates the slewing motion of the EVLA for use in 
determining optimized schedules. Based largely on the equivalent NRAO SSS code. 
"""

__author__ = "Seth Bruzewski"
__email__ = "bruzewskis@unm.edu"
__license__ = "GPL"

# Built in Modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord,Angle
from tqdm import tqdm, trange
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from SimEVLA import EquitorialToHorizontal as E2H, PointToPointTime as P2PT


def write_OPT(scheds, scs, priorities):
    '''
    Take in a list of schedules in schedmaster format and output appropriate
    VLA OPT schedule files
    '''
        
    # Templates for scan
    preamble1 = ('VERSION; 4;\n'
                 'SRC-CAT; 20A-439, VLA;\n'
                 'HDWR-CAT; NRAO Defaults, 20A-439;\n\n')
    
    preamble2 = ('SCHED-BLOCK; {schedBlockName}; {schedulingType}; '
                 '{iterationCount}; {date}; {timeOfDay}; {shadowLimit}; '
                 '{shadowCalcConfiguration}; {initTeleAz}; {initTeleEl}; '
                 '{avoidSunrise?}; {avoidSunset?}; {windApi}; '
                 '{commentsToOperator};\n\n')
    
    # Template for loops
    loopline = ('  LOOP-START; {loopName}; {loopIters}; {bracketed}; '
                '{comments};\n')
    
    loopend = '  LOOP-END;\n'
    
    # Template for actual scan lines
    scanline = ('  STD; {scanName}; {sourceName}; {resourceName}; {timeType}; '
                '{time}; {antennaWrap}; {applyRefPtg}; {applyPhase}; '
                '{recordOnMark5}; {allowOverTop}; {use10HzNoise}; '
                '{scanIntents}; {comments};\n')
    
    # Templates for sources
    sourceline = ('{sourceName}; {groupNames}; {coordinateSystem}; {epoch};'
                  '{lonCent}; {latCent}; {velocityRefFrame};'
                  '{velocityConvention}; {velocity}; {calibrator};\n')
    
    # Keep track of any sources we use
    sources = {} # name -> blocklist,RA,DEC
    
    # Loop on each block
    for i in trange(len(scheds)):
        
        # Calculate some things we want to know ahead of time
        # LST of first scan
        lst0 = scs[i].ra
        lst_str = lst0.to_string('h', pad=True)[:6]
        lst0 = Angle(lst0.hourangle//(5/60)*(5/60), unit=u.hourangle)
        
        # Figure out LST start range
        lst_min = lst0-60/60*u.hourangle
        lst_min_str = lst_min.wrap_at('360d').to_string('h', sep=':', pad=True)[:5]
        lst_max = lst0+60/60*u.hourangle
        lst_max_str = lst_max.wrap_at('360d').to_string('h', sep=':', pad=True)[:5]
        
        # open block file
        block_name = 'block_{0}_{2}{1:03d}'
        name_str = 'outputs/' + block_name + '.optScan'
        scanfile = name_str.format(lst_str, i+1, priorities[i])
        bf = open(scanfile, 'w')
        
        # Write preamble 1
        bf.write(preamble1)
        
        # Define preamble 2 params
        p2_parms = {'schedBlockName': block_name.format(lst_str, i+1, priorities[i]),
                    'schedulingType': 'Dynamic',
                    'iterationCount': 1,
                    'date': Time.now().iso[:-4]+', 2099-12-31 23:59:59',
                    'timeOfDay': lst_min_str+'-'+lst_max_str+',',
                    'shadowLimit': 0.0,
                    'shadowCalcConfiguration': 'C',
                    'initTeleAz': 225.0,
                    'initTeleEl': 35.0,
                    'avoidSunrise?': 'N',
                    'avoidSunset?': 'N',
                    'windApi': 'Any',
                    'commentsToOperator': ''}
        bf.write(preamble2.format(**p2_parms))
        
        # Convenient dictionary format
        sched = scheds[i]
        
        # Loop on each source
        for scan in sched:
            
            # Typical scan case
            if 'contents' not in scan:        
                # Write to line
                bf.write(scanline.format(**scan ))
                
                name = scan['sourceName']
                is_scan = scan['scanName'] == ''
                if is_scan and name not in sources:
                    sources[name] = {'ra': scan['ra'], 
                                     'dec': scan['dec']}
            else:
                
                bf.write(loopline.format(**scan))
                
                for lscan in scan['contents']:
                    bf.write('  '+scanline.format(**lscan))
                    
                    name = lscan['sourceName']
                    is_tgt = 'ObsTgt' in lscan['scanIntents']
                    if is_tgt and name not in sources:
                        sources[name] = {'ra': lscan['ra'],
                                         'dec': lscan['dec']}
                    
                bf.write(loopend)
                
        # close block file
        bf.close()
        
    
    # Open source list file
    sourcefile = 'outputs/EmptyFillers.pst'
    sf = open(sourcefile, 'w')
    for source in sorted(sources.keys()):
        # Define sourceline for source
        source_parms = {'sourceName': source,
                        'groupNames': '',
                        'coordinateSystem': 'Equatorial',
                        'epoch': 'J2000',
                        'lonCent': sources[source]['ra'],
                        'latCent': sources[source]['dec'],
                        'longRange': '',
                        'latRange': '',
                        'velocityRefFrame': '',
                        'velocityConvention': '',
                        'velocity': '',
                        'calibrator': ''}
        # Write
        sf.write(sourceline.format(**source_parms))
        
    # close source list file
    print('Source file written to', sourcefile)
    sf.close()
    
    return None

if __name__=='__main__':
    test = Table.read('testingblock.fits')
    testing = [test, test]
    test_ha = [0]*len(testing)
    write_OPT(testing, test_ha)
    
    