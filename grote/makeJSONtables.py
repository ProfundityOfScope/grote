#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 23:47:11 2021

@author: bruzewskis
"""

from dataclasses import dataclass, field
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import Union
import xmltodict
import json

if __name__ == '__main__':

    # Generate the source list in a less obnoxious way
    out = {'description': 'nraoDefaultSources', 'version': '1.1',
           'aliases': {}, 'gaincals': {}, 'sources': {}}
    
    # Set up gaincals
    for c in list('ABCD'):
        out['gaincals'][c] = {}
        for b in list('PLCXUKQ'):
            out['gaincals'][c][b] = {}
            for q in list('PSWCX'):
                out['gaincals'][c][b][q] = []
    
    with open('static/DefaultSources.xml', 'r') as fp:

        data = xmltodict.parse(fp.read(), process_namespaces=False)
        sources = data['sss:sourceCatalog']['sss:sources']['sss:source']

        cnt = 0
        for source in sources:

            # Find stuff we care about
            name = source['@name']
            aliases = source['sss:alias'] if 'sss:alias' in source else []
            aliasid = f'{cnt:04d}'
            out['sources'][aliasid] = {'name': name, 'notes': {}}

            # Write all names to the out
            out['aliases'][name] = aliasid
            for alias in aliases:
                out['aliases'][alias] = aliasid

            # Find center and write it
            subsources = source['sss:subsources']['sss:subsource']
            center = subsources['sss:simpleSkyPosition']
            lat = float(center['sss:latitude']['@value'])  # in sec
            lon = float(center['sss:longitude']['@value'])  # in arcsec
            sc = SkyCoord(lon/88400*360, lat/3600, unit='deg,deg')
            out['sources'][aliasid]['position'] = sc.to_string('hmsdms')

            # Write notes
            for udv in source['sss:userDefinedValues']['sss:userDefinedValue']:
                qual = udv['@value']
                out['sources'][aliasid]['notes'][udv['@key']] = qual
                
                if udv['@value'] in 'PSWCX':
                    config = udv['@key'][15]
                    band = udv['@key'][8]
                        
                    out['gaincals'][config][band][qual].append(aliasid)

            cnt += 1
        
        for c in list('ABCD'):
            for b in list('PLCXUKQ'):
                codes = out['gaincals'][c][b]
                codes['PS'] = codes['P'] + codes['S']
                codes['PSW'] = codes['PS'] + codes['W']
    

    # Generate the resource list in a less obnoxious way
    with open('static/DefaultResources.xml', 'r') as fp:

        data = xmltodict.parse(fp.read(), process_namespaces=False)

        out['resources'] = {}
        for entry in data['resourceCatalog']['entries']['resource']:

            name = entry['@name']
            out['resources'][name] = entry['notes']

    aliases = out['aliases']
    out['fluxcals'] = [aliases[c] for c in ['3C286', '3C48', '3C147', '3C138']]
    
    out['polcals'] = {'categoryA': [aliases[c] for c in ['3C48',       '3C138',      '3C286']],
                      'categoryB': [aliases[c] for c in ['J0359+5057', 'J0555+3948', 'J0854+2006', 
                                                         'J0927+3902', 'J1310+3220', 'J2136+0041', 
                                                         'J2202+4216','J2253+1608']],
                      'categoryC': [aliases[c] for c in ['3C84',      '3C147',       'J0713+4349',
                                                         'J1407+2827', 'J2355+4950']],
                      'categoryD': [aliases[c] for c in ['J0029+3456', 'J0111+3906', 'J0410+7656',
                                                         'J1035+5628', 'J1148+5924', 'J1400+6210',
                                                         'J1815+6127', 'J1823+7938', 'J1944+5448',
                                                         'J1945+7055', 'J2022+6136', 'J0022+0014',
                                                         'J0318+1628', 'J0329+2756', 'J1326+3154']]}


    json.dump(out, open('static/Defaults.json', 'w'))
