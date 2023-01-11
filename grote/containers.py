#!/usr/bin/env python3
"""
Containers script is the thing you're looking at.

There's more of a docstring here.
"""

import numpy as np
from astropy.coordinates import SkyCoord, Angle, EarthLocation
import astropy.units as u
from astropy.time import Time

from dataclasses import dataclass, field
from typing import Union
import json

import pkgutil

__author__ = 'Seth Bruzewski'
__copyright__ = "TBD"
__credits__ = ["Frank Schinzel", "Greg Taylor"]
__license__ = "TBD"
__version__ = "1.0.0"
__maintainer__ = "Seth Bruzewski"
__email__ = "bruzewskis@unm.edu"
__status__ = "Development"

_resources = json.loads(pkgutil.get_data('grote',
                                         'static/DefaultResources.json'))
_sources = json.loads(pkgutil.get_data('grote',
                                       'static/DefaultSources.json'))


@dataclass
class Resource:
    """
    Class theoretically describes a VLA resource.

    just to text.
    """

    name: str

    is_default: bool = False
    defaults = None

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

    @classmethod
    def from_name(cls, name):

        # Get relevant catalog
        resources = _resources['resources']

        # Search through defaults
        if name in resources:

            # Grab source and create it
            return cls(name, is_default=True)
        else:
            raise ValueError(f'{name} is not in default source catalog')


@dataclass
class Source:
    """Theoretically describes a source to be pointed at."""

    name: str
    coord: SkyCoord

    is_default: bool = False

    group_names: str = ''
    coord_system: str = 'Equatorial'
    epoch: str = 'J2000'
    velocity_frame: str = ''
    velocity_conv: str = ''
    velocity: str = ''
    calibrator: str = ''

    def __post_init__(self):
        self.lon = self.coord.ra.to_string('hour', sep=':', pad=True)
        self.lat = self.coord.dec.to_string(sep=':', pad=True, alwayssign=True)

    def __str__(self) -> str:
        pretty_pos = self.coord.to_string('hmsdms')
        return f'<Source: {self.name}, Position: ({pretty_pos})>'

    def __repr__(self) -> str:
        return self.__str__()

    @classmethod
    def from_name(cls, name):

        # Get relevant catalog
        sourcekeys = _sources['sourcekeys']

        # Search through defaults
        if name in sourcekeys:

            # Grab source and create it
            defsource = _sources['sources'][sourcekeys[name]]
            return cls(defsource['name'],
                       SkyCoord(defsource['position']),
                       is_default=True)
        else:
            raise ValueError(f'{name} is not in default source catalog')

    @classmethod
    def nearest_calibrator(cls, coord, config, band):

        qualstr = f'Quality.{band}-band.{config}-cfg'

        sourcemap = _sources['sources']
        calnotes = {k: sourcemap[k]['notes'] for k in sourcemap}
        goodcals = {k: v[qualstr] for k, v in calnotes.items() if qualstr in v}
        bettercals = [k for k in goodcals if goodcals[k] == 'P'
                      or goodcals[k] == 'S']

        # Get data from relevant catalog
        pos = SkyCoord([sourcemap[k]['position'] for k in bettercals])

        seps = coord.separation(pos).deg
        close = np.argmin(seps)
        cname = sourcemap[bettercals[close]]['name']

        return cls(cname,
                   pos[close],
                   is_default=True)

    def _write_str(self):

        # Templates for sources
        sourceline = ('{name}; {group_names}; {coord_system}; '
                      '{epoch}; {lon}; {lat}; {velocity_frame}; '
                      '{velocity_conv}; {velocity}; {calibrator};\n')

        # Quickly parse vars
        sourcevars = {}
        for k, v in vars(self).items():
            if isinstance(v, bool):
                sourcevars[k] = {True: 'Y', False: 'N'}[v]
            else:
                sourcevars[k] = v

        outstr = sourceline.format(**sourcevars)

        return outstr


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

    time_type: str = 'DUR'
    antenna_wrap: str = 'CW'
    apply_ref_ptg: bool = False
    apply_phase: bool = False
    record_on_mark5: bool = False
    allow_over_top: bool = False
    use_10Hz_noise: bool = True
    comments: str = ''

    def __post_init__(self):

        # Parse source string if we need to
        if isinstance(self.source, str):
            self.source = Source.from_name(self.source)

        # Parse resource string if we need to
        if isinstance(self.resource, str):
            self.resource = Resource.from_name(self.resource)

    def __str__(self):
        return self.pretty()

    def __repr__(self):
        pretty_pos = self.source.coord.to_string()
        return f'<Scan: {self.name}, Pos: ({pretty_pos})>'

    def pretty(self, depth=1, indent=0):
        pretty_time = self.time.to_string()
        return '\t'*indent + f'- {self.name} ({pretty_time})\n'

    def _write_str(self, known_sources):

        # Template for actual scan lines
        scanline = ('  STD; {name}; {source.name}; {resource.name}; '
                    '{time_type}; {time}; {antenna_wrap}; {apply_ref_ptg}; '
                    '{apply_phase}; {record_on_mark5}; {allow_over_top}; '
                    '{use_10Hz_noise}; {intents}; {comments};\n')

        # Quickly parse vars
        scanvars = {}
        for k, v in vars(self).items():
            if isinstance(v, bool):
                scanvars[k] = {True: 'Y', False: 'N'}[v]
            if isinstance(v, u.Quantity):
                # Check if its a time
                if v.unit.is_equivalent(u.hour):
                    time_as_angle = Angle(v.to(u.hour))
                    scanvars[k] = time_as_angle.to_string(sep=':', pad=True)
            else:
                scanvars[k] = v

        outstr = scanline.format(**scanvars)

        # Add any new sources to the table
        if self.source not in known_sources and not self.source.is_default:
            known_sources.append(self.source)

        return outstr


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

    is_bracketed: bool = False
    comments: str = ''

    def pretty(self, depth: int = 1, indent: int = 0) -> str:
        out = '\t'*indent + f'o {self.name} ({self.repeat}x)\n'
        for item in self.scans:
            if depth > 0:
                out += item.pretty(depth-1, indent+1)
        return out

    def _write_str(self, known_sources):

        # Template for loops
        loopline = ('  LOOP-START; {name}; {repeat}; {is_bracketed}; '
                    '{comments};\n')

        loopend = '  LOOP-END;\n'

        # Quickly parse vars
        loopvars = {}
        for k, v in vars(self).items():
            if isinstance(v, bool):
                loopvars[k] = {True: 'Y', False: 'N'}[v]
            else:
                loopvars[k] = v

        outstr = loopline.format(**loopvars)
        for scan in self.scans:
            outstr += scan._write_str(known_sources)
        outstr += loopend

        return outstr

@dataclass
class Block:
    """
    This class theoretically describes a block which contains loops or sources.
    """

    name: str
    scans: list[Union[Source, Loop, Mosaic]]
    config: str
    start_time: Union[Time, str]

    window: float = 1
    sched_type: str = 'Dynamic'
    iter_count: int = 1
    date: str = ''
    time_of_day: str = ''
    shadow_limit: int = 0
    init_tele_az: float = 225
    init_tele_el: float = 35
    avoid_sunrise: bool = False
    avoid_sunset: bool = False
    max_wind: float = 100
    max_phase: float = 175
    comments: str = ''

    def __post_init__(self):
        # Set the wind stuff
        self.wind_api = f'w={self.max_wind},p={self.max_phase}'

        # Convert start time to a time object
        self.loc = EarthLocation.of_site('VLA')
        if isinstance(self.start_time, str):
            self.start_time = Time(self.start_time, location=self.loc)

        # Set up date and time_of_day based on start time
        if self.sched_type == 'Fixed':
            # For fixed, we break up the provided iso time
            isotime = self.start_time.iso
            self.date = isotime[:10]
            self.time_of_day = isotime[11:19]
        elif self.sched_type == 'Dynamic':
            # Figure out start and end of window
            time1 = (self.start_time - self.window/2 * u.hour)
            time2 = (self.start_time + self.window/2 * u.hour)
            self.date = time1.iso

            # Set time of day to start and end LSTs
            lst1 = time1.sidereal_time('apparent').to_string(sep=':', pad=True)
            lst2 = time2.sidereal_time('apparent').to_string(sep=':', pad=True)
            self.time_of_day = f'{lst1[:5]}-{lst2[:5]}'

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

        # TODO: let this be able to handle loops and stuff
        ras = [s.source.coord.ra.wrap_at('180d').deg for s in self.scans]
        neworder = np.argsort(ras)

        self.scans = self.scans[neworder]

    def simulate(self) -> float:
        """Simulate project runtime."""
        time = 0 * u.s

        currentsc = SkyCoord('0h +60d')

        for i, item in enumerate(self.scans):
            if isinstance(item, Scan):

                # Slew to target
                slew = currentsc.separation(item.source.coord).deg
                slewtime = 0.5 * slew * u.s
                currentsc = item.source.coord

                # Add slew and on-source
                time += item.time + slewtime
            elif isinstance(item, Loop):
                for subscan in item.scans:

                    # Slew to target
                    slew = currentsc.separation(subscan.source.coord).deg
                    slewtime = 0.5 * slew * u.s
                    currentsc = subscan.source.coord

                    # Add slew and on-source
                    time += subscan.time + slewtime
            else:
                # TODO: implement mosaics
                pass

        return time

    def _write_str(self, known_sources) -> str:

        preamble2 = ('SCHED-BLOCK; {name}; {sched_type}; '
                     '{iter_count}; {date}; {time_of_day}; {shadow_limit}; '
                     '{config}; {init_tele_az}; {init_tele_el}; '
                     '{avoid_sunrise}; {avoid_sunset}; {wind_api}; '
                     '{comments};\n\n')

        # Quickly parse vars
        blockvars = {}
        for k, v in vars(self).items():
            if isinstance(v, bool):
                blockvars[k] = {True: 'Y', False: 'N'}[v]
            else:
                blockvars[k] = v

        # Loop through contents
        outstr = preamble2.format(**blockvars)
        for item in self.scans:
            outstr += item._write_str(known_sources)

        return outstr

    def write(self, blockname, srcname, overwrite=True) -> None:

        # If called directly, write out a block file
        preamble1 = ('VERSION; 4;\n'
                     'SRC-CAT; EXAMPLE, VLA;\n'
                     'HDWR-CAT; NRAO Defaults, EXAMPLE;\n\n')

        known_sources = []
        with open(blockname, 'w') as bp:
            bp.write(preamble1)
            bp.write(self._write_str(known_sources))

        # If called directly, we need to also write source file
        sourcestr = ''
        for source in known_sources:
            sourcestr += source._write_str()

        with open(srcname, 'w') as sp:
            for source in known_sources:
                sp.write(source._write_str())


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

    b1 = Block('Block1', [s1, s2, l1], 'C', '2023-01-10 15:59:43')
    b1 = Block('Block2', [s1, l2, s2], 'C', '2023-01-10 20:59:13')

    p1 = Project('Project1', [b1, b1])
    print(p1.pretty(3))

    print(p1.simulate())
    print(b1.simulate())

    return p1


if __name__ == '__main__':
    test_project = make_test_project()
