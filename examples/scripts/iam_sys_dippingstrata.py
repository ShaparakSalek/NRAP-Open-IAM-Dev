"""
Example demonstrating the use of the DippingStratigraphy component. Here, 
the component is used to calculate unit thicknesses and depths at one 
location.

Example of run:
$ python iam_sys_dippingstrata.py
"""

import sys
import os
import logging
import numpy as np
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, DippingStratigraphy


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)
    
    # The reference location at which the thickness parameters apply
    locXRef = 1200
    locYRef = 1200
    
    # The location where thicknesses and depths will be calculated
    locX = 3000
    locY = 3000
    
    numberOfShaleLayers = 3
    
    # Thicknesses at the reference location
    shale_thicknesses = [750, 950, 200]
    aquifer_thicknesses = [200, 200]
    reservoir_thickness = 150
    
    strike = 315
    dip = 5
    dipDirection = 'NE'
    
    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    
    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)
    
    # Add stratigraphy component
    strata = sm.add_component_model_object(DippingStratigraphy(
        name='strata', parent=sm, locXRef=locXRef, locYRef=locYRef,
        dipDirection=dipDirection, locX=locX, locY=locY))
    
    strata.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    strata.add_par('reservoirThickness', value=reservoir_thickness, vary=False)
    
    for shaleRef in range(numberOfShaleLayers):
        strata.add_par('shale{}Thickness'.format(shaleRef + 1), 
                       value=shale_thicknesses[shaleRef], vary=False)
        
        if (shaleRef + 1) < numberOfShaleLayers:
            strata.add_par('aquifer{}Thickness'.format(shaleRef + 1), 
                           value=aquifer_thicknesses[shaleRef], vary=False)
    
    strata.add_par('strike', value=strike, vary=False)
    strata.add_par('dip', value=dip, vary=False)
    
    # Only use get_thickness_obs_names() and get_depth_obs_names() after the 
    # numberOfShaleLayers parameter has been assigned.
    thickness_obs = strata.get_thickness_obs_names()
    
    # The thicness and depths observations are only produced in the first time 
    # step, so use index=[0] when adding these observations.
    for ob_nm in thickness_obs:
        strata.add_obs(ob_nm, index=[0])
    
    depth_obs = strata.get_depth_obs_names()
    
    for ob_nm in depth_obs:
        strata.add_obs(ob_nm, index=[0])
    
    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    # Collect thickness observations
    reservoirThickness = sm.collect_observations_as_time_series(
        strata, 'reservoirThickness', indices=[0])
    shale1Thickness = sm.collect_observations_as_time_series(
        strata, 'shale1Thickness', indices=[0])
    aquifer1Thickness = sm.collect_observations_as_time_series(
        strata, 'aquifer1Thickness', indices=[0])
    shale2Thickness = sm.collect_observations_as_time_series(
        strata, 'shale2Thickness', indices=[0])
    aquifer2Thickness = sm.collect_observations_as_time_series(
        strata, 'aquifer2Thickness', indices=[0])
    shale3Thickness = sm.collect_observations_as_time_series(
        strata, 'shale3Thickness', indices=[0])
    
    print('Unit thicknesses (m):')
    print('    reservoirThickness: ', reservoirThickness)
    print('    shale1Thickness: ', shale1Thickness)
    print('    aquiferThickness: ', aquifer1Thickness)
    print('    shale2Thickness: ', shale2Thickness)
    print('    aquifer2Thickness: ', aquifer2Thickness)
    print('    shale3Thickness: ', shale3Thickness)
    print('')
    
    # Collect depth observations. Because these are only produced during the 
    # first timestep, use indices=[0]. Otherwise, it will look for an observation 
    # named 'strata.reservoirDepth_#', where '_#' is the time step number. The 
    # DippingStratigraphy component does not produce observations with  
    # names that include the time step number.
    reservoirDepth = sm.collect_observations_as_time_series(
        strata, 'reservoirDepth', indices=[0])
    shale1Depth = sm.collect_observations_as_time_series(
        strata, 'shale1Depth', indices=[0])
    aquifer1Depth = sm.collect_observations_as_time_series(
        strata, 'aquifer1Depth', indices=[0])
    shale2Depth = sm.collect_observations_as_time_series(
        strata, 'shale2Depth', indices=[0])
    aquifer2Depth = sm.collect_observations_as_time_series(
        strata, 'aquifer2Depth', indices=[0])
    shale3Depth = sm.collect_observations_as_time_series(
        strata, 'shale3Depth', indices=[0])
    
    print('Depths (m) to the base of each unit:')
    print('    reservoirDepth: ', reservoirDepth)     
    print('    shale1Depth: ', shale1Depth)
    print('    aquifer1Depth: ', aquifer1Depth)
    print('    shale2Depth: ', shale2Depth)
    print('    aquifer2Depth: ', aquifer2Depth)
    print('    shale3Depth: ', shale3Depth)
    
    reservoirMidDepth = sm.collect_observations_as_time_series(
        strata, 'reservoirMidDepth', indices=[0])
    shale1MidDepth = sm.collect_observations_as_time_series(
        strata, 'shale1MidDepth', indices=[0])
    aquifer1MidDepth = sm.collect_observations_as_time_series(
        strata, 'aquifer1MidDepth', indices=[0])
    shale2MidDepth = sm.collect_observations_as_time_series(
        strata, 'shale2MidDepth', indices=[0])
    aquifer2MidDepth = sm.collect_observations_as_time_series(
        strata, 'aquifer2MidDepth', indices=[0])
    shale3MidDepth = sm.collect_observations_as_time_series(
        strata, 'shale3MidDepth', indices=[0])
    
    print('Depths (m) to the middle of each unit:')
    print('    reservoirMidDepth: ', reservoirMidDepth)     
    print('    shale1MidDepth: ', shale1MidDepth)
    print('    aquifer1MidDepth: ', aquifer1MidDepth)
    print('    shale2MidDepth: ', shale2MidDepth)
    print('    aquifer2MidDepth: ', aquifer2MidDepth)
    print('    shale3MidDepth: ', shale3MidDepth)
    
    reservoirTopDepth = sm.collect_observations_as_time_series(
        strata, 'reservoirTopDepth', indices=[0])
    shale1TopDepth = sm.collect_observations_as_time_series(
        strata, 'shale1TopDepth', indices=[0])
    aquifer1TopDepth = sm.collect_observations_as_time_series(
        strata, 'aquifer1TopDepth', indices=[0])
    shale2TopDepth = sm.collect_observations_as_time_series(
        strata, 'shale2TopDepth', indices=[0])
    aquifer2TopDepth = sm.collect_observations_as_time_series(
        strata, 'aquifer2TopDepth', indices=[0])
    shale3TopDepth = sm.collect_observations_as_time_series(
        strata, 'shale3TopDepth', indices=[0])
    
    print('Depths (m) to the top of each unit:')
    print('    reservoirTopDepth: ', reservoirTopDepth)     
    print('    shale1TopDepth: ', shale1TopDepth)
    print('    aquifer1TopDepth: ', aquifer1TopDepth)
    print('    shale2TopDepth: ', shale2TopDepth)
    print('    aquifer2TopDepth: ', aquifer2TopDepth)
    print('    shale3TopDepth: ', shale3TopDepth)
    
