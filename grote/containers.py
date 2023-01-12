#!/usr/bin/env python3
"""
Containers script is the thing you're looking at.

There's more of a docstring here.
"""

from dataclasses import dataclass, field
from typing import Union
import json
import pkgutil
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time

from grote.SimEVLA import MotionSimulator

# Some things these functions will regularly access
_ppath = Path(__file__).parent
_resources = json.load(open(_ppath / 'static/Defaults.json', 'r'))


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

        # Search through defaults
        if name in _resources['aliases']:
            alias = _resources['aliases'][name]

            # Grab source and create it
            source = _resources['sources'][alias]
            fname = source['name']
            fpos = SkyCoord(source['position'])
            return cls(fname, fpos, is_default=True)
        else:
            raise ValueError(f'{name} is not in default source catalog')

    @classmethod
    def nearest_calibrator(cls, coord, config, band, quality='PS'):

        band = 'C' if band == 'S' else band
        cals = _resources['gaincals'][config][band][quality]

        # Get data from relevant catalog
        pos = SkyCoord([_resources['sources'][a]['position'] for a in cals])

        seps = coord.separation(pos).deg
        closest = np.argmin(seps)

        cname = _resources['sources'][cals[closest]]['name']
        cpos = pos[closest]

        return cls(cname, cpos, is_default=True)

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
    """Scan container class."""

    name: str
    source: Union[Source, str]
    resource: Union[Resource, str]
    time: u.Quantity
    intents: str

    time_type: str = 'DUR'
    antenna_wrap: str = 'No Preference'
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
        return '\t'*indent + f'- {self.name} ({pretty_time}, {self.intents})\n'

    def _write_str(self, known_sources, slew_time):

        # Template for actual scan lines
        scanline = ('STD; {name}; {source.name}; {resource.name}; '
                    '{time_type}; {duration}; {antenna_wrap}; '
                    '{apply_ref_ptg}; {apply_phase}; {record_on_mark5}; '
                    '{allow_over_top}; {use_10Hz_noise}; {intents}; '
                    '{comments};\n')

        # Quickly parse vars
        dur = slew_time + self.time
        dur_ang = Angle(dur.to(u.hour)).to_string(sep=':', pad=True)
        scanvars = {'duration': dur_ang}
        for k, v in vars(self).items():
            if isinstance(v, bool):
                scanvars[k] = 'Y' if v else 'N'
            else:
                scanvars[k] = v

        outstr = scanline.format(**scanvars)

        # Add any new sources to the table
        if self.source not in known_sources and not self.source.is_default:
            known_sources.append(self.source)

        return outstr


@dataclass
class Mosaic:
    """Mosaic container class."""

    name: str


@dataclass
class Loop:
    """Loop container class."""

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

    def _write_str(self, known_sources, slew_times):

        # Template for loops
        loopline = ('LOOP-START; {name}; {repeat}; {is_bracketed}; '
                    '{comments};\n')

        loopend = '  LOOP-END;\n'

        # Quickly parse vars
        loopvars = {}
        for k, v in vars(self).items():
            if isinstance(v, bool):
                loopvars[k] = 'Y' if v else 'N'
            else:
                loopvars[k] = v

        outstr = loopline.format(**loopvars)
        for scan, slew in zip(self.scans, slew_times):
            outstr += '\t' + scan._write_str(known_sources, slew)
        outstr += loopend

        return outstr


@dataclass
class Block:
    """Block container class."""

    name: str
    scans: list[Union[Source, Loop, Mosaic]]
    config: str
    start_time: Union[Time, str]

    window: float = 1
    valid_schedule = False
    slews = None

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
            lst1 = time1.sidereal_time('apparent').hms
            lst1_str = f'{lst1.h:02.0f}:{5*(lst1.m//5):02.0f}'

            time2 = (self.start_time + self.window/2 * u.hour)
            lst2 = time2.sidereal_time('apparent').hms
            lst2_str = f'{lst2.h:02.0f}:{5*(lst2.m//5):02.0f}'

            self.date = time1.iso[:10]
            self.time_of_day = f'{lst1_str}-{lst2_str}'

    def __str__(self) -> str:
        return self.pretty()

    def __repr__(self) -> str:
        pretty_scans = [s.name for s in self.scans]
        return f'<Block: {self.name}, Scans: {pretty_scans}>'

    def pretty(self, depth: int = 1, indent: int = 0) -> str:
        """Generate a pretty recursive depiction."""
        out = '\t'*indent + f'> {self.name}\n'
        for item in self.scans:
            if depth > 0:
                out += item.pretty(depth-1, indent+1)
        return out

    def optimize(self, threshold: float = 0.01) -> None:

        # TODO: let this be able to handle loops and stuff
        pass

    def simulate(self, plot=False) -> float:
        """Simulate project runtime."""

        currenttime = self.start_time
        loc = EarthLocation.of_site('VLA')
        sim = MotionSimulator()

        self.slews = []
        times = []
        els = []
        intents = []
        for i, item in enumerate(self.scans):
            if isinstance(item, Scan):

                # Slew to  target
                frame = AltAz(obstime=currenttime, location=loc)
                tgt = item.source.coord.transform_to(frame)
                slewtime = sim.moveTo(tgt.az, tgt.alt)

                # Record
                self.slews.append(slewtime)
                times.append(currenttime)
                els.append(tgt.alt.deg)
                intents.append(item.intents)
                currenttime += item.time + slewtime

                # Update after following
                frame = AltAz(obstime=currenttime, location=loc)
                tgt = item.source.coord.transform_to(frame)
                _ = sim.moveTo(tgt.az, tgt.alt)

            elif isinstance(item, Loop):
                loopslews = []
                for _ in range(item.repeat):
                    for subscan in item.scans:

                        # Slew to  target
                        frame = AltAz(obstime=currenttime, location=loc)
                        tgt = subscan.source.coord.transform_to(frame)
                        slewtime = sim.moveTo(tgt.az, tgt.alt)

                        # Record
                        loopslews.append(slewtime)
                        times.append(currenttime)
                        els.append(tgt.alt.deg)
                        intents.append(subscan.intents)
                        currenttime += subscan.time + slewtime

                        # Update after following
                        frame = AltAz(obstime=currenttime, location=loc)
                        tgt = subscan.source.coord.transform_to(frame)
                        _ = sim.moveTo(tgt.az, tgt.alt)

                if item.is_bracketed:
                    # Slew to target
                    scan0 = item.scans[0]

                    # Slew to  target
                    frame = AltAz(obstime=currenttime, location=loc)
                    tgt = scan0.source.coord.transform_to(frame)
                    slewtime = sim.moveTo(tgt.az, tgt.alt)

                    # Record
                    loopslews.append(slewtime)
                    times.append(currenttime)
                    els.append(tgt.alt.deg)
                    intents.append(scan0.intents)
                    currenttime += scan0.time + slewtime

                    # Update after following
                    frame = AltAz(obstime=currenttime, location=loc)
                    tgt = scan0.source.coord.transform_to(frame)
                    _ = sim.moveTo(tgt.az, tgt.alt)

                # Log the slews for the loop
                self.slews.append(loopslews)

            else:
                # TODO: implement mosaics
                pass

        if np.min(els) > 15:
            self.valid_schedule = True

        if plot:
            import matplotlib.pyplot as plt

            cmap = {}
            for i, un in enumerate(np.unique(intents)):
                cmap[un] = f'C{i}'
                plt.scatter(0, -10, c=f'C{i}', label=un)

            for i in range(len(times)):
                runtime = (times[i].mjd - times[0].mjd) * 24
                plt.scatter(runtime, els[i], c=cmap[intents[i]])
            plt.ylim(0, 90)
            plt.legend()
            plt.grid()
            plt.title(f'Schedule Starting at {self.time_of_day[:5]} LST')
            plt.show()

        # Calculate and return total time
        total_time = currenttime - self.start_time
        return total_time

    def _write_str(self, known_sources) -> str:

        # Check if the schedule has been validated
        if not self.valid_schedule:
            totaltime = self.simulate(plot=True).to('hour')
            print('Schedule has been validated, totals to:', totaltime)

        preamble2 = ('SCHED-BLOCK; {name}; {sched_type}; '
                     '{iter_count}; {date}; {time_of_day}; {shadow_limit}; '
                     '{config}; {init_tele_az}; {init_tele_el}; '
                     '{avoid_sunrise}; {avoid_sunset}; {wind_api}; '
                     '{comments};\n\n')

        # Quickly parse vars
        blockvars = {}
        for k, v in vars(self).items():
            if isinstance(v, bool):
                blockvars[k] = 'Y' if v else 'N'
            else:
                blockvars[k] = v

        # Loop through contents
        outstr = preamble2.format(**blockvars)
        for item, slew in zip(self.scans, self.slews):
            outstr += item._write_str(known_sources, slew)

        return outstr

    def write(self, blockname, srcname, overwrite=True) -> None:
        """Public write function for people."""
        # If called directly, write out a block file
        preamble1 = ('VERSION; 4;\n'
                     'SRC-CAT; EXAMPLE, VLA;\n'
                     'HDWR-CAT; NRAO Defaults;\n\n')

        known_sources = []
        with open(blockname, 'w') as bp:
            bp.write(preamble1)
            bp.write(self._write_str(known_sources))

        # If called directly, we need to also write source file
        with open(srcname, 'w') as sp:
            for source in known_sources:
                sp.write(source._write_str())


@dataclass
class Project:
    """Project container class."""

    name: str = 'Default'
    blocks: list[Block] = field(default_factory=list)

    def __str__(self) -> str:
        """Magic string method."""
        return self.pretty()

    def __repr__(self) -> str:
        """Magic repr method."""
        pretty_blocks = [b.name for b in self.blocks]
        return f'<Project: {self.name}, Blocks: {pretty_blocks}>'

    def pretty(self, depth: int = 1, indent: int = 0) -> str:
        """Generate a pretty recursive depiction."""
        out = '\t'*indent + self.name + '\n'
        for item in self.blocks:
            if depth > 0:
                out += item.pretty(depth-1, indent+1)
        return out

    def simulate(self) -> float:
        """Simulate project runtime."""
        time = 0 * u.s
        for block in self.blocks:
            blocktime = block.simulate().to(u.hour)
            time += blocktime

            print(block.name, 'total runtime:', blocktime)

        return time

    def write(self, schedname: str, srcname: str, overwrite: bool = False) -> None:

        # If called directly, write out a block file
        preamble1 = ('VERSION; 4;\n'
                     'SRC-CAT; {name}, VLA;\n'
                     'HDWR-CAT; NRAO Defaults;\n\n')

        # Write out big file
        known_sources = []
        for i,block in enumerate(self.blocks):
            with open(schedname.replace('.', f'_b{i}.'), 'w') as bp:
                bp.write(preamble1.format(**vars(self)))
                bp.write(block._write_str(known_sources))

        # If called directly, we need to also write source file
        with open(srcname, 'w') as sp:
            for source in known_sources:
                sp.write(source._write_str())

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

    s1 = Scan('Scan1', '3C286', xr, 4*u.min, 'setAtnGain')
    s2 = Scan('Scan1', '3C286', xr, 1*u.min, 'CalBP')
    s3 = Scan('Scan1', '3C286', xr, 3*u.min, 'CalPolLeak')
    s4 = Scan('Science', Source('Source2', SkyCoord('02h24m32s', '88d15\'30"')),
              lr, 10*u.min, 'ObsTgt')

    c1 = Source.nearest_calibrator(s4.source.coord, 'A', 'L', 'PS')
    s5 = Scan('Calib', c1, lr, 1*u.min, 'CalGain')

    l1 = Loop('Loop1', 5, [s5, s4], is_bracketed=True)

    b1 = Block('Block1', [s1, s2, s3, l1], 'C', '2023-01-10 08:59:43')

    p1 = Project('Project1', [b1])
    print(p1.pretty(3))
    dummy = b1.simulate(True)

    return p1


if __name__ == '__main__':
    test_project = make_test_project()
