#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
schedmaster.py: This code identifies candidates for followup with 12 hours of
further VLA time on FERMI sources, then helps to create a schedule, testing
that it's possible. If everything looks good, we send the schedule to the
to be converted to OPT format and write it to disk. Note that this code serves
also as a prototype for the pythonic OPT project. 
"""

__author__ = "Seth Bruzewski"
__email__ = "bruzewskis@unm.edu"
__license__ = "GPL"

# Third Party Modules
import numpy as np
from astropy.coordinates import SkyCoord,Angle
import astropy.units as u

# Custom Modules
from SimEVLA import EquitorialToHorizontal as E2H
from SimEVLA import PointToPointTime as P2P
from SimEVLA import MotionSimulator as vlams

def find_nearest_cal(sc, calibrators, config):
    '''
    For an input coordinate, returns the nearest high quality VLA calibrator
    for the specified band and array configuration. 
    
    Args
    =======
    sc (SkyCoord) - The sky coordinate to search around for a calibrator.
    band (string) - The observational band.
    config (string) - The observational array configuration.
    
    Returns
    =======
    best_sc (SkyCoord) - The skycoordinate of the nearest calibrator.
    name (string) - The J2000 epoch name for the returned calibrator.
    
    Raises
    =======
    None
    '''
    
    # Narrow down to good quality
    config_qual = calibrators[config+'_qual']
    good_qual = np.logical_or(config_qual=='P', config_qual=='S')
    cals = calibrators[good_qual]
    
    # Only grab things that won't cause a flip
    zenith = 34.07875
    cal_decs = Angle(cals['dec'])
    if sc.dec.deg < zenith:
        cals = cals[cal_decs.deg < zenith]
    else:
        cals = cals[cal_decs.deg > zenith]
    
    # Convert to SkyCoords
    sc_cals = SkyCoord(cals['ra'], cals['dec'])
    
    # Find minimum distance calibrator
    seps = sc.separation(sc_cals).deg
    near = np.argmin(seps)
    best_sc = sc_cals[near]
    name = cals['name'][near]
    
    return best_sc, name

def worst_setup_slew(az, el):
    '''
    This function finds the worst possible setup slew, assuming the telescope
    begins the observation at the worst possible wrap end. 
    
    Args
    =======
    az (Angle) - The target azimuth of the first source
    el (Angle) - The target elevation of the first source
    
    Returns
    =======
    worst_time (Quantity) - The time the worst possible slew will take
    
    Raises
    =======
    None
    '''
    
    # HARD CODING SLEW SPEEDS
    vel_az = 40/60 * u.Unit('deg/s')
    vel_el = 20/60 * u.Unit('deg/s')
    
    # HARD CODING WRAP LIMITS
    az_ends = np.array([-85, 445])
    el_ends = np.array([0,90])
    
    # Figure out max in each direction
    az_dist = max(abs(az_ends-az.deg)) * u.deg
    el_dist = max(abs(el_ends-el.deg)) * u.deg
    
    time = max([ az_dist/vel_az, el_dist/vel_el])
    
    return time

def to_dur_fmt(t):
    '''
    Convert to format duration likes, input is time object
    '''
    time_hour = t.to('hour').value
    time_str = Angle(time_hour, unit=u.hourangle).to_string('h', pad=True)
    
    return time_str

def build_line(scanName='', sourceName='', resourceName='L16f5DC-realfast',
               time=2*u.min, antennaWrap='', scanIntents='ObsTgt, ', 
               ra=Angle('0d'), dec=Angle('0d')):
    
    line = {'scanName': scanName,
            'sourceName': sourceName,
            'resourceName': resourceName,
            'timeType': 'DUR',
            'time': to_dur_fmt(time),
            'antennaWrap': antennaWrap,
            'applyRefPtg': 'N',
            'applyPhase': 'N',
            'recordOnMark5': 'N',
            'allowOverTop': 'N',
            'use10HzNoise': 'Y',
            'scanIntents': scanIntents,
            'comments': '',
            'ra': ra.to_string('h', sep=':', pad=True),
            'dec': dec.to_string('deg', sep=':', pad=True, alwayssign=True)}
    
    return line

def schedule_block(source_sc, source_name, calibrators, verbose=True):
    '''
    For input points and start LST, run a simulated observation to see if 
    the particular schedule is possible. This will produce some figures
    displaying the hour angle and elevation of the telescope over time, 
    highlighting in particular any areas where the telescope would go beyond
    certian limits. The function returns the total time the observation would
    take
    
    Args
    =======
    block (Table) - Target list in the FERMI format. The needed columns are 
        essentially just RAJ2000, DEJ2000, and Source_Name.
    mosaics (dict) - The mosaic dictionary maping each source name to its
        mosaic coordinates
    source_scan (float/Quantity) - The amount of time to spend on each 
        pointing. Default value is 30 and default units are seconds.
    cal_scan (float/Quantity) - The amount of time to spend on phase
        calibrators. Default value is 60 and default units are seconds.
    ha_offset (Quantity) - The hour angle of the first source when the
        observation will begin. For instance, if RA[0]=1h and ha_offset=-1h,
        then the observation will start at LST=0h. Defaults to 0h, such that
        things start when the first source is on the meridian. Default units
        are hourangle
    verbose (bool) - Whether to generate HA and El plots of the projected
        schedule and produce explanatory text.
    band (string) - The band which the observation will be performed in. This
        is used to search for calibrators.
    config (string) - The array configuration when the observation will be
        performed. This is used to search for calibrators.
    
    Returns
    =======
    sched_table (Table) - A table summarizing the observation. If the user
        likes how it went, this table can be passed to OPTSCHEDULER for
        formatting.
    
    Raises
    =======
    None
    '''
    
    ############### SET UP STUFF ###############
    # Check units
    loop_iterations = 5
    source_scan = 12 * u.min
    cal_scan = 90 * u.s
    config = 'C'
    
    # Convenience function for time->ang
    time_to_ang = lambda t : Angle(t.to('s').value/3600, unit=u.hourangle)
    
    ############### CALIBRATOR LISTS ###############
    # Flux calibrators
    flux_cal_names = ['0521+166=3C138', '0542+498=3C147', '1331+305=3C286',
                      '0137+331=3C48']
    flux_cals = SkyCoord(['05h21m09.886021s +16d38\'22.051220"', # 3C138
                          '05h42m36.137916s +49d51\'07.233560"', # 3C147
                          '13h31m08.287984s +30d30\'32.958850"', # 3C286
                          '01h37m41.299431s +33d09\'35.132990"'])# 3C48
    
    # Polarization angle calibrators
    ang_cal_names = ['J0359+5057', 'J2253+1608']
    ang_cals = SkyCoord(['03h59m +50d57\'', '22h53m +16d08\''])
    
    
    # Polarization leakage calibrators
    leak_cal_names = ['J0319+4130', 'J0542+4951'	, 'J0713+4349', 
                      'J1407+2827', 'J2355+4950']
    leak_cals = SkyCoord(['03h19m +41d30\'', # J0319+4130
                          '05h42m +49d51\'', # J0542+4951 = 3C147
                          '07h13m +43d49\'', # J0713+4349
                          '14h07m +28d27\'', # J1407+2827
                          '23h55m +49d50\''])# J2355+4950
    
    ############### INITIALIZE TIMEKEEPING ###############
    # Figure out when to start
    lst = source_sc.ra.wrap_at('360d')
    
    ############### DETERMINE BEST FLUX CALIBRATOR ###############
    flux_az, flux_el = E2H(flux_cals.ra, flux_cals.dec, lst)
    best_ind = np.argmax(flux_el)
    
    # Values for best flux cal
    best_flux_name = flux_cal_names[best_ind]
    best_flux_sc = flux_cals[best_ind]
    
    # What is the worst possible slew we might have to do
    worst_slew = worst_setup_slew(flux_az[best_ind], flux_el[best_ind])
    setup_slew = 5*u.min
    setup_atten = 1*u.min
    setup_req = 30*u.s
    scan_tgt1 = 1*u.min
    scan_tgt2 = 1*u.min
    scan_polr = 2*u.min
    
    ############### SETUP SCANS ###############
    # Define scheduler keeper
    sched = []
    
    # Slew scan
    start_lst = lst
    end_lst = start_lst + time_to_ang(setup_slew)
    slew_line = build_line(scanName='slew', sourceName=best_flux_name, 
                           resourceName='X band pointing', time=setup_slew,
                           antennaWrap='', scanIntents='SetAtnGain, ', 
                           ra=best_flux_sc.ra, dec=best_flux_sc.dec)
    sched.append(slew_line)
    
    # Atten scan
    start_lst = end_lst
    end_lst = start_lst + time_to_ang(setup_atten)
    atten_line = build_line(scanName='atten', sourceName=best_flux_name, 
                            time=setup_atten, antennaWrap='', 
                            scanIntents='SetAtnGain, ', ra=best_flux_sc.ra, 
                            dec=best_flux_sc.dec)
    sched.append(atten_line)
    
    # Req scan
    start_lst = end_lst
    end_lst = start_lst + time_to_ang(setup_req)
    req_line = build_line(scanName='req', sourceName=best_flux_name, 
                          time=setup_req, antennaWrap='', 
                          scanIntents='SetAtnGain, ', ra=best_flux_sc.ra, 
                          dec=best_flux_sc.dec)
    sched.append(req_line)
    
    ############### CALIBRATION SCANS ###############
    # Flux target scan 1
    start_lst = end_lst
    end_lst = start_lst + time_to_ang(scan_tgt1)
    flux_time = scan_tgt1 + worst_slew - setup_slew - setup_atten - setup_req
    flux_tgt1_line = build_line(sourceName=best_flux_name, time=flux_time, 
                                scanIntents='Calgain, ', ra=best_flux_sc.ra, 
                                dec=best_flux_sc.dec)
    sched.append(flux_tgt1_line)
    
    # Apply different flux and pol calibration strategies
    if best_flux_name=='0542+498=3C147':
        # Scan 3c147 as a leakage calibrator
        start_lst = end_lst
        end_lst = start_lst + time_to_ang(scan_tgt2)
        flux_tgt2_line = build_line(sourceName=best_flux_name, time=scan_tgt2, 
                                    scanIntents='CalBP, CalFlux, CalPolLeak, ', 
                                    ra=best_flux_sc.ra, dec=best_flux_sc.dec)
        sched.append(flux_tgt2_line)
        
        # Find best polarization angle calibrator
        ang_az, ang_el = E2H(ang_cals.ra, ang_cals.dec, lst)
        best_ind = np.argmax(ang_el)
        best_ang_name = ang_cal_names[best_ind]
        best_ang_sc = ang_cals[best_ind]
        
        # Scan pol ang calibrator
        slew_to_ang = P2P(best_flux_sc, best_ang_sc, end_lst)
        start_lst = end_lst
        end_lst = start_lst + time_to_ang(scan_polr + slew_to_ang)
        ang_time = scan_polr + slew_to_ang
        ang_tgt_line = build_line(sourceName=best_ang_name, time=ang_time, 
                                  scanIntents='CalGain, CalPolAng, ', 
                                  ra=best_ang_sc.ra, dec=best_ang_sc.dec)
        sched.append(ang_tgt_line)
        
        # Record final position
        end_cal_sc = best_ang_sc
        
    elif best_flux_name=='0137+331=3C48':
        # Scan flux cal without any polarization info
        start_lst = end_lst
        end_lst = start_lst + time_to_ang(scan_tgt2)
        flux_tgt2_line = build_line(sourceName=best_flux_name, time=scan_tgt2, 
                                    scanIntents='CalBP, CalFlux, ', 
                                    ra=best_flux_sc.ra, dec=best_flux_sc.dec)
        sched.append(flux_tgt2_line)
        
        # Find best polarization angle calibrator
        ang_az, ang_el = E2H(ang_cals.ra, ang_cals.dec, lst)
        best_ind = np.argmax(ang_el)
        best_ang_name = ang_cal_names[best_ind]
        best_ang_sc = ang_cals[best_ind]
        
        # Scan pol ang calibrator
        slew_to_ang = P2P(best_flux_sc, best_ang_sc, end_lst)
        start_lst = end_lst
        end_lst = start_lst + time_to_ang(scan_polr + slew_to_ang)
        ang_time = scan_polr + slew_to_ang
        ang_tgt_line = build_line(sourceName=best_ang_name, time=ang_time, 
                                  scanIntents='CalGain, CalPolAng, ', 
                                  ra=best_ang_sc.ra, dec=best_ang_sc.dec)
        sched.append(ang_tgt_line)
        
        # Find best polarization leakage calibrator
        leak_az, leak_el = E2H(leak_cals.ra, leak_cals.dec, lst)
        best_ind = np.argmax(leak_el)
        best_leak_name = leak_cal_names[best_ind]
        best_leak_sc = leak_cals[best_ind]
        
        # Scan pol leak calibrator
        slew_to_leak = P2P(best_ang_sc, best_leak_sc, end_lst)
        start_lst = end_lst
        end_lst = start_lst + time_to_ang(scan_polr + slew_to_leak)
        leak_time = scan_polr + slew_to_leak
        leak_tgt_line = build_line(sourceName=best_leak_name, time=leak_time, 
                                   scanIntents='CalGain, CalPolLeak, ', 
                                   ra=best_leak_sc.ra, dec=best_leak_sc.dec)
        sched.append(leak_tgt_line)
        
        # Record final position
        end_cal_sc = best_leak_sc
        
    else:
        # Scan flux cal as angle calibrator
        start_lst = end_lst
        end_lst = start_lst + time_to_ang(scan_tgt2)
        flux_tgt2_line = build_line(sourceName=best_flux_name, time=scan_tgt2, 
                                    scanIntents='CalBP, CalFlux, CalPolAng, ', 
                                    ra=best_flux_sc.ra, dec=best_flux_sc.dec)
        sched.append(flux_tgt2_line)
        
        # Find best polarization leakage calibrator
        leak_az, leak_el = E2H(leak_cals.ra, leak_cals.dec, lst)
        best_ind = np.argmax(leak_el)
        best_leak_name = leak_cal_names[best_ind]
        best_leak_sc = leak_cals[best_ind]
        
        # Scan pol leak calibrator
        slew_to_leak = P2P(best_flux_sc, best_leak_sc, end_lst)
        start_lst = end_lst
        end_lst = start_lst + time_to_ang(scan_polr + slew_to_leak)
        leak_time = scan_polr + slew_to_leak
        leak_tgt_line = build_line(sourceName=best_leak_name, time=leak_time, 
                                   scanIntents='CalGain, CalPolLeak, ', 
                                   ra=best_leak_sc.ra, dec=best_leak_sc.dec)
        sched.append(leak_tgt_line)
        
        # Record final position
        end_cal_sc = best_leak_sc
    
    # Gain calibrator entry (before-loop)
    cal_sc, cal_name = find_nearest_cal(source_sc, calibrators, config)
    slew_to_cal = P2P(end_cal_sc, cal_sc, end_lst)
    start_lst = end_lst
    end_lst = start_lst + time_to_ang(cal_scan + slew_to_cal)
    cal_time = cal_scan + slew_to_cal
    cali_tgt_line = build_line(sourceName=cal_name, time=cal_time, 
                               scanIntents='CalGain, ', 
                               ra=cal_sc.ra, dec=cal_sc.dec)
    sched.append(cali_tgt_line)
        
    ############### SCAN LOOP ENTRIES ###############
    # Tell it a loop is coming
    loop = {'loopName': 'Targets', 
            'loopIters': loop_iterations, 
            'bracketed': 'N',
            'comments': '',
            'contents': []}
    
    # Calculate theoretical time we should spend on source
    loops = 5
    cal_to_src = P2P(cal_sc, source_sc, end_lst) # Assume its the same
    flux_time = (end_lst - lst).hour * u.hour
    non_source_time = (2*cal_to_src + cal_scan)*loops
    opt_source_time = (1.5*u.hour - flux_time - non_source_time)/loops
    
    source_scan = max([source_scan, opt_source_time])
    
    # Target scan
    source_time = source_scan + cal_to_src
    source_line = build_line(sourceName=source_name, time=source_time, 
                             scanIntents='ObsTgt, ', 
                             ra=source_sc.ra, dec=source_sc.dec)
    loop['contents'].append(source_line)

    # Calibrator scan
    cal_time = cal_scan + cal_to_src
    cal_line = build_line(sourceName=cal_name, time=cal_time, 
                          scanIntents='CalGain, ', 
                          ra=cal_sc.ra, dec=cal_sc.dec)
    loop['contents'].append(cal_line)
    
    # Add loop to scans
    sched.append(loop)
    
    return sched
    
if __name__=='__main__':
    
    print('Hello!')