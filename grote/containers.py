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

@dataclass
class Resource:
    '''
    This class theoretically describes a VLA resource. 
    '''
    name: str
    library: str = field(init=False)
    
    def __post_init__(self):
        default_resources = {'Lx32'}
        if self.name in default_resources:
            self.library = 'NRAO_Default'
        else:
            self.library = 'MyResources'

@dataclass
class Source:
    '''
    This class theoretically describes a source to be pointed at
    '''
    
    name: str
    coord: SkyCoord
    resource: Resource
    dur: u.Quantity = 30*u.s
    
    def __str__(self):
        return self.pretty()
    
    def __repr__(self):
        pretty_pos = self.coord.to_string()
        return f'<Scan: {self.name}, Pos: ({pretty_pos})>'
    
    def pretty(self, depth=1, indent=0):
        pretty_time = self.dur.to_string()
        return '\t'*indent + f'- {self.name} ({pretty_time})\n'
    
    def separation(self, other):
        return self.coord.separation(other.coord).deg
    
@dataclass
class Mosaic:
    '''
    This class theoretically describes a source to be pointed at
    '''
    
    name: str
    coord: SkyCoord
    resource: Resource
    dur: u.Quantity = 30*u.s
    num_pointings : int = 7 # this will probably end up post init
    
    def __str__(self):
        return self.pretty()
    
    def __repr__(self):
        pretty_pos = self.coord.to_string()
        return f'<Scan: {self.name}, Pos: ({pretty_pos})>'
    
    def pretty(self, depth=1, indent=0):
        return '\t'*indent + f'= {self.name} [{self.num_pointings}p]\n'
    
    def separation(self, other):
        return self.coord.separation(other.coord).deg
    
@dataclass
class Loop:
    '''
    This class theoretically describes a loop which contains sources
    '''
    name: str
    repeat: int
    scans: list[Source, Mosaic]
    
    def __str__(self):
        return self.pretty()
    
    def __repr__(self):
        pretty_scans = [ s.name for s in self.scans ]
        return f'<Loop: {self.name}, Scans: {pretty_scans}>'
    
    def pretty(self, depth=1, indent=0):
        out = '\t'*indent + f'o {self.name} [x{self.repeat}]\n'
        for item in self.scans:
            if depth>0:
                out += item.pretty(depth-1, indent+1)
        return out

@dataclass
class Block:
    '''
    This class theoretically describes a block which contains loops or sources
    '''
    name: str
    start_time: str
    scans: list[Union[Source, Loop, Mosaic]]
    
    def __str__(self):
        return self.pretty()
    
    def __repr__(self):
        pretty_scans = [ s.name for s in self.scans ]
        return f'<Block: {self.name}, Scans: {pretty_scans}>'
    
    def pretty(self, depth : int = 1, indent : int = 0) -> str:
        out = '\t'*indent + f'> {self.name}\n'
        for item in self.scans:
            if depth>0:
                out += item.pretty(depth-1, indent+1)
        return out
    
    @classmethod
    def from_file(cls, filename : str, stype : str = 'fits'):
        '''
        This method should read a block from a file. It should be flexible 
        enough to handle a csv or fits table with the right columns 
        (Name,RA,DEC,Intent,etc...) or just a file right from the OPT. Can 
        either have the user entry fits/csv/opt or we can guess it
        '''
        return cls()
    
@dataclass
class Project:
    '''
    This class theoretically describes a project which contains blocks
    '''
    name: str = 'Default'
    blocks: list[Block] = field(default_factory=list)
    
    def __str__(self):
        return self.pretty()
    
    def __repr__(self):
        pretty_blocks = [ b.name for b in self.blocks ]
        return f'<Project: {self.name}, Blocks: {pretty_blocks}>'
    
    def pretty(self, depth=1, indent=0):
        out = '\t'*indent + self.name + '\n'
        for item in self.blocks:
            if depth>0:
                out += item.pretty(depth-1, indent+1)
        return out
    
    @classmethod
    def from_xml(cls, filename : str):
        '''
        Not implemented yet, will eventually return a full constructed project
        which one can then edit as they like
        '''
        return cls()
    
    def write(self, filename : str, style : str = 'xml', 
              clobber : bool = False) -> bool:
        '''
        Not implemented yet, will eventually write out the file either as XML
        or as all the relevant text files one would need
        '''
        return True
    
    def simulate(self) -> float:
        '''
        Simple implementation of timing, assuming no slew time, just adds up
        the time
        '''
        time = 0
        for block in self.blocks:
            for scan in block.scans:
                if isinstance(scan, Loop):
                    loop_time = 0
                    for sub_scan in scan.scans:
                        loop_time += sub_scan.dur
                    time += scan.repeat * loop_time
                else:
                    time += scan.dur
                    
        return time
    
def make_test_project():
    Lx32 = Resource('Lx32')
        
    s1 = Source('Source1', SkyCoord('01h01m01s','01d01\'01"'), Lx32, 5*u.min)
    s2 = Source('Source2', SkyCoord('02h02m02s','02d02\'02"'), Lx32)
    s3 = Source('Source3', SkyCoord('03h03m03s','03d03\'03"'), Lx32)
    m1 = Mosaic('Mosaic1', SkyCoord('03h03m03s','03d03\'03"'), Lx32, 7)
    l1 = Loop('Loop1', 2, [s2,s3,m1])
    b1 = Block('Block1', '00:00:00.00', [s1,l1])
    
    s4 = Source('Source4', SkyCoord('01h01m01s','01d01\'01"'), Lx32, 5*u.min)
    s5 = Source('Source5', SkyCoord('02h02m02s','02d02\'02"'), Lx32)
    s6 = Source('Source6', SkyCoord('03h03m03s','03d03\'03"'), Lx32)
    m2 = Mosaic('Mosaic2', SkyCoord('03h03m03s','03d03\'03"'), Lx32, 7)
    m3 = Mosaic('Mosaic3', SkyCoord('03h03m03s','03d03\'03"'), Lx32, 7)
    l2 = Loop('Loop2', 2, [s3,s5,s6, m3])
    b2 = Block('Block2', '00:00:00.00', [s4, m2, l2])
    
    
    p1 = Project('Project1', [b1,b2])
    print(p1.pretty(3))
    # print(p1)
    # print(p1.simulate())
    return p1

if __name__=='__main__':
    test_project = make_test_project()