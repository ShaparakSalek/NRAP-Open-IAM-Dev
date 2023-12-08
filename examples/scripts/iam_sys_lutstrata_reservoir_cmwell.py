"""
This example demonstrates the use of LookupTableStratigraphy component. This 
component is used to produce unit thicknesses at one location through spatial 
interpolation of the data stored in a .csv file. At the location used, the 
unit thicknesses produced are linked to AnalyticalReservoir and CementedWellbore 
components.

When using a CementedWellbore component with LookupTableStratigraphy, the 
depthRatio parameter needs to be set with the add_par_linked_to_composite_obs() 
method. This method is similar to the add_composite_par() method, in that the 
user provides an expression that will be evaluated to solve for the value of a 
parameter. While add_composite_par() calculates a parameter value based on other 
parameter values, the add_par_linked_to_composite_obs() method can calculate a 
parameter value based on a combination of parameters and observations.

Example of run:
$ python iam_sys_lutstrata_reservoir_cmwell.py
"""

import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, StratigraphyDataInterpolator,
                     LookupTableStratigraphy, AnalyticalReservoir, 
                     CementedWellbore)


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'examples', 'Control_Files', 
                                  'input_data', 'ex38a'])
    
    data_file = 'stratigraphy.csv'

    # Define keyword arguments of the system model
    num_years = 25
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    
    numberOfShaleLayers = 3
    
    # In the .csv file used, the center of the points is at x = 0 km, y = 0 km
    injX = 0
    injY = 0
    
    # The index for the main aquifer receicing the leakage from the cemented wellbore
    AquiferIndex = 2
    
    # The index for the thief zone aquifer, attenuating the leakage from the 
    # cemented wellbore to the main aquifer
    ThiefZoneIndex = 1
    
    # Using x and y values outside of the range covered in the .csv file will 
    # cause LookupTableStratigraphy component to fail at the spatial interpolation.
    locX = 1767.77
    locY = 1767.77
    
    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)
    
    intpr = sm.add_interpolator(StratigraphyDataInterpolator(
        name='int1', parent=sm,
        header_file_dir=file_directory,
        data_file=data_file), 
        intr_family='stratigraphy')

    # Add stratigraphy component
    luts = sm.add_component_model_object(LookupTableStratigraphy(
        name='luts', parent=sm, intr_family='stratigraphy', locX=locX, locY=locY))
    
    luts.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    
    thickness_obs = luts.get_thickness_obs_names()
    
    for ob_nm in thickness_obs:
        luts.add_obs(ob_nm, index=[0])
        luts.add_obs_to_be_linked(ob_nm)
    
    depth_obs = luts.get_depth_obs_names()
    
    for ob_nm in depth_obs:
        luts.add_obs(ob_nm, index=[0])
        luts.add_obs_to_be_linked(ob_nm)
    
    ares = sm.add_component_model_object(AnalyticalReservoir(
        name='ares', parent=sm, injX=injX, injY=injY, locX=locX, locY=locY))
    
    ares.add_par('logResPerm', value=-13.0, vary=False)
    ares.add_par('injRate', value=1.0, vary=False)
    # This value should not be too low relative to the locX and locY used
    ares.add_par('reservoirRadius', value=10000, vary=False)
    ares.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    
    # Link the unit thickness parameters to the LookupTableStratigraphy component
    for ob_nm in thickness_obs:
        ares.add_par_linked_to_obs(ob_nm, luts.linkobs[ob_nm])
    
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')
    
    cmw = sm.add_component_model_object(CementedWellbore(
        name='cmw', parent=sm))
    
    # Add parameters of cemented wellbore component
    cmw.add_par('wellRadius', value=0.025, vary=False)
    cmw.add_par('logWellPerm', value=-13.0, vary=False)
    cmw.add_par('logThiefPerm', value=-12.25, vary=False)
    
    cmw.add_par_linked_to_obs('wellDepth', luts.linkobs['shale1Depth'])
    
    # The parameter depthRatio is a function of unit thicknesses. When using a 
    # LookupTableStratigraphy component, unit thicknesses are prouced as observations. 
    # Set depthRatio with add_par_linked_to_composite_obs().
    depth_ratio = '(luts.aquifer{tz}Thickness/2'.format(tz=ThiefZoneIndex)
    for il in range(ThiefZoneIndex + 1, numberOfShaleLayers):
        depth_ratio += ' + luts.shale{il}Thickness + '.format(
            il=il) + 'luts.aquifer{il}Thickness'.format(il=il)
        
    depth_ratio += ' + luts.shale{nsl}Thickness)'.format(
        nsl=numberOfShaleLayers)
    
    depth_ratio += '/cmw.wellDepth'
    
    cmw.add_par_linked_to_composite_obs('depthRatio', depth_ratio)
    
    cmw.add_composite_par('initPressure',
        expr='ares.datumPressure + cmw.wellDepth*cmw.g*ares.brineDensity')
    
    cmw.add_obs('brine_aquifer1')
    cmw.add_obs('brine_aquifer2')
    cmw.add_obs('brine_atm')
    cmw.add_obs('CO2_aquifer1')
    cmw.add_obs('CO2_aquifer2')
    cmw.add_obs('CO2_atm')
    cmw.add_obs('mass_CO2_aquifer1')
    cmw.add_obs('mass_CO2_aquifer2')
    
    cmw.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    cmw.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])
    
    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    
    # Collect unit depth observations
    reservoirDepth = sm.collect_observations_as_time_series(luts, 'reservoirDepth', 
                                                            indices=[0])
    shale1Depth = sm.collect_observations_as_time_series(luts, 'shale1Depth', 
                                                         indices=[0])
    shale2Depth = sm.collect_observations_as_time_series(luts, 'shale2Depth', 
                                                         indices=[0])
    shale3Depth = sm.collect_observations_as_time_series(luts, 'shale3Depth', 
                                                         indices=[0])
    aquifer1Depth = sm.collect_observations_as_time_series(luts, 'aquifer1Depth', 
                                                           indices=[0])
    aquifer2Depth = sm.collect_observations_as_time_series(luts, 'aquifer2Depth', 
                                                           indices=[0])
    
    print('Depths (m) to the bottom of each unit:')
    print('    reservoirDepth: ', reservoirDepth)     
    print('    shale1Depth: ', shale1Depth)
    print('    aquiferDepth: ', aquifer1Depth)
    print('    shale2Depth: ', shale2Depth)
    print('    aquifer2Depth: ', aquifer2Depth)
    print('    shale3Depth: ', shale3Depth)
    
    # Make the figures
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': 10}
    plt.rc('font', **font)
    
    pressure = sm.collect_observations_as_time_series(ares, 'pressure')
    CO2saturation = sm.collect_observations_as_time_series(ares, 'CO2saturation')
    
    time_array_yrs = time_array / 365.25
    
    plt.figure(1, figsize=(12, 8))
    
    ax = plt.subplot(2,1,1)
    ax.plot(time_array_yrs, pressure, color='k', linewidth=2)
    ax.set_ylabel('Pressure (Pa)', fontsize=12, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
    
    ax = plt.subplot(2,1,2)
    ax.plot(time_array_yrs, CO2saturation, color='r', linewidth=2)
    ax.set_ylabel('CO$_2$ Saturation', fontsize=12, fontweight='bold')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
    
    CO2_aquifer1 = sm.collect_observations_as_time_series(cmw, 'CO2_aquifer1')
    CO2_aquifer2 = sm.collect_observations_as_time_series(cmw, 'CO2_aquifer2')
    CO2_atm = sm.collect_observations_as_time_series(cmw, 'CO2_atm')
    
    plt.figure(2, figsize=(9, 7))
    
    ax=plt.gca()
    ax.plot(time_array_yrs, CO2_aquifer1, color='b', linewidth=2, 
            label='To Aquifer 1')
    ax.plot(time_array_yrs, CO2_aquifer2, color='r', linestyle='--', 
            linewidth=2.5, label='To Aquifer 2')
    ax.plot(time_array_yrs, CO2_atm, color='k', linestyle=':', 
            linewidth=2.5, label='To Atmosphere')
    ax.set_ylabel('CO$_2$ Leakage Rate to Aquifer (kg/s)', fontsize=12, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
    plt.legend(fontsize=12)
    
    mass_CO2_aquifer1 = sm.collect_observations_as_time_series(cmw, 
                                                               'mass_CO2_aquifer1')
    mass_CO2_aquifer2 = sm.collect_observations_as_time_series(cmw, 
                                                               'mass_CO2_aquifer2')
    
    plt.figure(3, figsize=(9, 7))
    
    ax=plt.gca()
    ax.plot(time_array_yrs, mass_CO2_aquifer1, color='b', linewidth=2, 
            label='Aquifer 1')
    ax.plot(time_array_yrs, mass_CO2_aquifer2, color='r', linestyle='--', 
            linewidth=2.5, label='Aquifer 2')
    ax.set_ylabel('CO$_2$ Leaked Mass (kg)', fontsize=12, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
    plt.legend(fontsize=12)
    
