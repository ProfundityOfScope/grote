#!/usr/bin/env python3
"""
Containers script is the thing you're looking at.

There's more of a docstring here.
"""

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

from dataclasses import dataclass, field
from typing import Union
import json

from importlib.resources import files, Package
import static

__author__ = 'Seth Bruzewski'
__copyright__ = "TBD"
__credits__ = ["Frank Schinzel", "Greg Taylor"]
__license__ = "TBD"
__version__ = "1.0.0"
__maintainer__ = "Seth Bruzewski"
__email__ = "bruzewskis@unm.edu"
__status__ = "Development"


@dataclass
class Resource:
    """
    Class theoretically describes a VLA resource.

    just to text.
    """

    name: str

    isdefault: bool = False

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.__str__()
    
    @classmethod
    def from_name(cls, name):

        # Open default catalog
        # test = files('static').join()
        with files('static').joinpath('DefaultResources.json').open('rb') as fp:
            dat = json.load(fp)
            resources = dat['resources']

            # Search through defaults
            if name in resources:

                # Grab source and create it
                return cls(name, isdefault=True)
            else:
                raise ValueError(f'{name} is not in default source catalog')


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
    
    @classmethod
    def nearest_calibrator(cls, coord, config, band):

        # Open default catalog
        with open('../static/DefaultSources.json', 'r') as fp:
            dat = json.load(fp)
            
            keys = list(dat['sources'].keys())
            pos = SkyCoord([ dat['sources'][k]['position'] for k in keys ])
            
            seps = coord.separation(pos).deg
            close = np.argmin(seps)
            
            return cls(dat['sources'][keys[close]],
                       pos[close],
                       isdefault=True)


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
        return '\t'*indent + f'- {self.name}\n'


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
    
    isbracked: bool = False
    comments: str = ''
    
    def pretty(self, depth: int = 1, indent: int = 0) -> str:
        out = '\t'*indent + f'o {self.name} ({self.repeat}x)\n'
        for item in self.scans:
            if depth > 0:
                out += item.pretty(depth-1, indent+1)
        return out


@dataclass
class Block:
    """
    This class theoretically describes a block which contains loops or sources
    """

    name: str
    start_time: str
    scans: list[Union[Source, Loop, Mosaic]]
    config: str

    lst: str = '00:00'
    schedtype: str = 'Dynamic'
    itercount: int = 1
    shadowlimit: int = 15
    initTeleAz: float = 225
    initTeleEl: float = 35
    avoidSunrise: bool = False
    avoidSunset: bool = False
    maxWind: float = 100
    maxPhase: float = 175

    def __str__(self) -> str:
        return self.pretty()

    def __repr__(self) -> str:
        pretty_scans = [s.name for s in self.scans]
        return f'<Block: {self.name}, Scans: {pretty_scans}>'

    def pretty(self, depth: int = 1, indent: int = 0) -> str:
        out = '\t'*indent + f'> {self.name}\n'
        for item in self.scans:
            if depth > 0:
                out += item.pretty(depth-1, indent+1)
        return out
    
    def optimize(self, threshold: float = 0.01) -> None:
        
        ras = [ s.source.coord.ra.wrap_at('180d').deg for s in self.scans ]
        neworder = np.argsort(ras)
        
        self.scans = self.scans[neworder]

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

    def write(self, filename, overwrite=True) -> None:
        pass


@dataclass
class Project:
    """
    This class theoretically describes a project which contains blocks
    """
    
    name: str = 'Default'
    blocks: list[Block] = field(default_factory=list)

    def __str__(self) -> str:
        return self.pretty()

    def __repr__(self) -> str:
        pretty_blocks = [b.name for b in self.blocks]
        return f'<Project: {self.name}, Blocks: {pretty_blocks}>'

    def pretty(self, depth: int = 1, indent: int = 0) -> str:
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

    def write(self, naming: str, overwrite: bool = False) -> None:
        for block in self.blocks:
            block.write(f'{naming}_{block.name}.blah')


def make_test_project():
    """
    Literally just a little testing function.

    Returns
    -------
    p1 : Project
        The generated project object.

    """
    xr = Resource.from_name('X-point')
    lr = Resource.from_name('L16f2A')

    s1 = Scan('Scan1', 'J0638+5933',
              xr, 15*u.min, 'target')
    s2 = Scan('Scan2', Source('Source2', SkyCoord('02h24m32s', '10d15\'30"')),
              lr, 15*u.min, 'target')
    s3 = Scan('Scan3', Source('Source3', SkyCoord('03h03m03s', '03d03\'03"')),
              lr, 20*u.min, 'target')

    c1 = Source.nearest_calibrator(s3.source.coord, 'A', 'L')
    s4 = Scan('CalScan', c1, lr, 5*u.min, 'cal')

    l1 = Loop('Loop1', 10, [s3, s4])
    l2 = Loop('Loop2', 5, [s4, s3])

    b1 = Block('Block1', '00:00:00.00', [s1, s2, l1], 'C')
    b1 = Block('Block2', '00:00:00.00', [s1, l2, s2], 'C')

    p1 = Project('Project1', [b1, b1])
    print(p1.pretty(3))

    print(p1.simulate())
    print(b1.simulate())

    return p1


if __name__ == '__main__':
    test_project = make_test_project()
