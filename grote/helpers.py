#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 00:55:10 2023

@author: bruzewskis
"""
import matplotlib.pyplot as plt
from matplotlib.cm import viridis
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, Angle, AltAz
from astropy.time import Time
import astropy.units as u
from time import time
import json
from pathlib import Path

# Some things these functions will regularly access
_ppath = Path(__file__).parent
_resources = json.load(open(_ppath / 'static/Defaults.json', 'r'))

def show_cals(starttime):
    loc = EarthLocation.of_site('VLA')
    time = Time(starttime, location=loc)

    data = _resources
    
    plt.figure(figsize=(10,10), dpi=192)
    
    plt.subplot2grid((2,2),(0,0), colspan=2, projection='mollweide')
    fluxcals = [data['sources'][aid]['position'] for aid in data['fluxcals']]
    fn = [data['sources'][aid]['name'] for aid in data['fluxcals']]
    fsc = SkyCoord(fluxcals, obstime=time, location=loc)
    fx = fsc.altaz.az.wrap_at('180d').rad
    fy = fsc.altaz.alt.rad
    for i in range(len(fsc)):
        plt.text(fx[i], fy[i], fn[i], c='C4', va='center', ha='center')
    plt.scatter(2*np.pi, 0, c='C4', label='Flux Cal')
    plt.grid()
    plt.legend()
    lst = time.sidereal_time('apparent').to_string(pad=True, sep=':', precision=0)
    plt.title(f'LST: {lst}')
    
    plt.subplot2grid((2,2),(1,0), colspan=2, projection='mollweide')
    for i,let in enumerate(list('DCBA')):
        cat = f'category{let}'
        polcals = [data['sources'][aid]['position'] for aid in data['polcals'][cat]]
        pn = [data['sources'][aid]['name'] for aid in data['polcals'][cat]]
        psc = SkyCoord(polcals, obstime=time, location=loc)
        px = psc.altaz.az.wrap_at('180d').rad
        py = psc.altaz.alt.rad
        for j in range(len(psc)):
            plt.text(px[j], py[j], pn[j], c=f'C{i}', va='center', ha='center')
    plt.scatter(2*np.pi, 0, c='C0', label='Secondary Leak Cals')
    plt.scatter(2*np.pi, 0, c='C1', label='Primary Leak Cals')
    plt.scatter(2*np.pi, 0, c='C2', label='Secondary Pol Cals')
    plt.scatter(2*np.pi, 0, c='C3', label='Primary Pol Cals')
    plt.grid()
    plt.legend(loc='lower right')
    plt.show()
    
def forecast(scs, starttime):
    loc = EarthLocation.of_site('VLA')
    time = Time(starttime, location=loc)

    data = _resources
    
    startlst = time.sidereal_time('apparent')
    diffs = np.linspace(0,1,100) * u.day
    
    times = time + diffs
    lsts = startlst.hour + np.linspace(0,24,100)
    
    frames = AltAz(obstime=times, location=loc)
    for i in range(len(scs)):
        aa = scs[i].transform_to(frames)
        plt.plot(lsts, aa.alt.deg, c=viridis(i/len(scs)))
    
    plt.grid()
    plt.ylim(0,90)
    plt.show()
    
if __name__=='__main__':
    sc = SkyCoord(np.random.normal(170,5,10),
                  np.random.normal(50, 5, 10), unit='deg,deg')
    forecast(sc, Time.now().iso)