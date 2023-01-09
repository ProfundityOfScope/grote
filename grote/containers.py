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
import json

@dataclass
class Resource:
    """
    Class theoretically describes a VLA resource.

    just to text.
    """

    name: str
    
    isdefault: bool = False
    
    @classmethod
    def from_name(cls, name):
        # Search through defaults
        if True:
            return cls(name, isdefault=True)
        else:
            raise ValueError(f'{name} is not in default resource catalog')


@dataclass
class Source:
    """Theoretically describes a source to be pointed at."""

    name: str
    coord: SkyCoord

    isdefault: bool = False

    groupNames: str = ''
    coordinateSystem: str = 'Equatorial'
    epoch: str = 'J2000'
    longRange: str = ''
    latRange: str = ''
    velocityRefFrame: str = ''
    velocityConvention: str = ''
    velocity: str = ''
    calibrator: str = ''

    def __str__(self) -> str:
        pretty_pos = self.coord.to_string('hmsdms')
        return f'<Source: {self.name}, Position: ({pretty_pos})>'

    def __repr__(self) -> str:
        return self.__str__()

    @classmethod
    def from_name(cls, name):

        # Open default catalog
        with open('../static/DefaultSources.json', 'r') as fp:
            dat = json.load(fp)
            sourcekeys = dat['sourcekeys']

            # Search through defaults
            if name in sourcekeys:

                # Grab source and create it
                defsource = dat['sources'][sourcekeys[name]]
                return cls(defsource['name'],
                           SkyCoord(defsource['position']),
                           isdefault=True)
            else:
                raise ValueError(f'{name} is not in default source catalog')


@dataclass
class Scan:
    """
    This is a scan, separate from a source
    """

    name: str
    source: Union[Source, str]
    resource: Union[Resource, str]
    time: u.Quantity
    intents: str

    timetype: str = 'DUR'
    antennaWrap: str = 'CW'
    applyRefPtg: bool = False
    applyPhase: bool = False
    recordOnMark5: bool = False
    allowOverTop: bool = False
    use10HzNoise: bool = True
    comments: str = ''

    def __post_init__(self):
        
        # Parse source string if we need to
        if isinstance(self.source, str):
            self.source = Source.from_name(self.source)
        
        # Parse resource string if we need to
        if isinstance(self.resource, str):
            self.resource = Resource.from_name(self.source)

    def __str__(self):
        return self.pretty()

    def __repr__(self):
        pretty_pos = self.coord.to_string()
        return f'<Scan: {self.name}, Pos: ({pretty_pos})>'

    def pretty(self, depth=1, indent=0):
        pretty_time = self.time.to_string()
        return '\t'*indent + f'- {self.name} ({self.source.name}) [{pretty_time}]\n'


@dataclass
class Mosaic:
    """
    This class theoretically describes a source to be pointed at
    """

    name: str


@dataclass
class Loop:
    """
    This class theoretically describes a loop which contains sources
    """

    name: str
    repeat: int
    scans: list[Scan, Mosaic]


@dataclass
class Block:
    """
    This class theoretically describes a block which contains loops or sources
    """

    name: str
    start_time: str
    scans: list[Union[Source, Loop, Mosaic]]

    def __str__(self):
        return self.pretty()

    def __repr__(self):
        pretty_scans = [s.name for s in self.scans]
        return f'<Block: {self.name}, Scans: {pretty_scans}>'

    def pretty(self, depth: int = 1, indent: int = 0) -> str:
        out = '\t'*indent + f'> {self.name}\n'
        for item in self.scans:
            if depth > 0:
                out += item.pretty(depth-1, indent+1)
        return out

    def simulate(self) -> float:
        """Simulate project runtime.

        A simple implementation of timing, assuming no slew time, just adds up
        the time.

        Returns
        -------
        time : float
            The time.
        """
        time = 0 * u.s
        for i, scan in enumerate(self.scans):
            if isinstance(scan, Scan):
                time += scan.time
            else:
                pass  # TODO: add loop and mosaic

        return time


@dataclass
class Project:
    """
    This class theoretically describes a project which contains blocks
    """
    
    name: str = 'Default'
    blocks: list[Block] = field(default_factory=list)

    def __str__(self):
        return self.pretty()

    def __repr__(self):
        pretty_blocks = [b.name for b in self.blocks]
        return f'<Project: {self.name}, Blocks: {pretty_blocks}>'

    def pretty(self, depth=1, indent=0):
        out = '\t'*indent + self.name + '\n'
        for item in self.blocks:
            if depth > 0:
                out += item.pretty(depth-1, indent+1)
        return out

    @classmethod
    def from_xml(cls, filename: str):
        """
        Not implemented yet, will eventually return a full constructed project
        which one can then edit as they like
        """
        return cls()

    def write(self, filename: str, style: str = 'xml',
              overwrite: bool = False) -> bool:
        """
        Not implemented yet, will eventually write out the file either as XML
        or as all the relevant text files one would need
        """
        return True

    def simulate(self) -> float:
        """Simulate project runtime.

        A simple implementation of timing, assuming no slew time, just adds up
        the time.

        Returns
        -------
        time : float
            The time.
        """
        time = 0 * u.s
        for block in self.blocks:
            time += block.simulate()

        return time


def make_test_project():
    """
    Literally just a little testing function.

    Returns
    -------
    p1 : Project
        The generated project object.

    """
    Lx32 = Resource.from_name('Lx32')

    s1 = Scan('Scan1', 'J0638+5933',
              Lx32, 15*u.min, 'target')
    s2 = Scan('Scan2', Source('Source2', SkyCoord('02h02m02s', '02d02\'02"')),
              Lx32, 15*u.min, 'target')
    s3 = Scan('Scan3', Source('Source3', SkyCoord('03h03m03s', '03d03\'03"')),
              Lx32, 20*u.min, 'target')

    b1 = Block('Block1', '00:00:00.00', [s1, s2, s3])

    p1 = Project('Project1', [b1, b1])
    print(p1.pretty(3))
    print(p1.simulate())
    print(b1.simulate())

    return p1


if __name__ == '__main__':
    test_project = make_test_project()