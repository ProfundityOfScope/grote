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
    out = {'description': 'nraoDefaultSources', 'version': '1.0.0',
           'sourcekeys': {}, 'sources': {}}
    with open('../static/DefaultSources.xml', 'r') as fp:
        
        data = xmltodict.parse(fp.read(), process_namespaces=False)
        sources = data['sss:sourceCatalog']['sss:sources']['sss:source']
        
        cnt = 0
        for source in sources:
            
            # Find stuff we care about
            name = source['@name']
            aliases = source['sss:alias'] if 'sss:alias' in source else []
            out['sources'][cnt] = {'name': name, 'notes': {}}
            
            # Write all names to the out
            out['sourcekeys'][name] = str(cnt)
            for alias in aliases:
                out['sourcekeys'][alias] = str(cnt)
                
            # Find center and write it
            subsources = source['sss:subsources']['sss:subsource']
            center = subsources['sss:simpleSkyPosition']
            lat = float(center['sss:latitude']['@value']) # in sec
            lon = float(center['sss:longitude']['@value']) # in arcsec
            sc = SkyCoord(lon/88400*360, lat/3600, unit='deg,deg')
            out['sources'][cnt]['position'] = sc.to_string('hmsdms')
            
            # Write notes
            for udv in source['sss:userDefinedValues']['sss:userDefinedValue']:
                out['sources'][cnt]['notes'][udv['@key']] = udv['@value']
                
            cnt += 1
            
        # Dump it in a nicer format
        json.dump(out, open('../static/DefaultSources.json', 'w'))
        
    # Generate the resource list in a less obnoxious way
    with open('../static/DefaultResources.xml', 'r') as fp:
        
        data = xmltodict.parse(fp.read(), process_namespaces=False)
        
        out = {'description': 'nraoDefaultResources', 'version': '1.0.0', 
               'resources': {}}
        for entry in data['resourceCatalog']['entries']['resource']:
            
            name = entry['@name']
            out['resources'][name] = {}
            out['resources'][name]['notes'] = entry['notes']
            
            for udv in entry['userDefinedValues']['userDefinedValue']:
                out['resources'][name]['notes'][udv['@key']] = udv['@value']
        
        # Dump it out
        json.dump(out, open('../static/DefaultResources.json', 'w'))