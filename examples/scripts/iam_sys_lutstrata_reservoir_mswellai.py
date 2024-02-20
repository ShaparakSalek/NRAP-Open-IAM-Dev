"""
This example demonstrates the use of LookupTableStratigraphy component. This 
component is used to produce unit thicknesses at one location through spatial 
interpolation of the data stored in a .csv file. At the location used, the 
unit thicknesses produced are linked to AnalyticalReservoir and MultisegmentedWellboreAI 
components. For example, the AnalyticalReservoir and MultisegmentedWellboreAI 
components have parameters like reservoirThickness, aquifer1Thickness, and 
shale3Thickness that are linked to the observations of the LookupTableStratigraphy 
component.

To run this example, the files for the Multisegmented Wellbore AI component must be 
downloaded. The files are not included in the download of NRAP-Open-IAM because they are 
large (about 7.2 GB). The two files can be downloaded from these links:
    http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/42882a99-a5b0-4d68-8f4c-cad4812358c4/download
    http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/8b786e28-0e52-40d4-84b5-23c25dc2d338/download

The zipped files must then be unzipped, with the resulting folders placed into this directory:
    src/openiam/components/models/wellbore/multisegmented_ai
    
Example of run:
$ python iam_sys_lutstrata_reservoir_mswellai.py
"""

import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam.components.iam_base_classes import SystemModel
from openiam.components.stratigraphy_data_interpolator import StratigraphyDataInterpolator
from openiam.components.lookup_table_stratigraphy_component import LookupTableStratigraphy
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.multisegmented_wellbore_ai_component import MultisegmentedWellboreAI

if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'examples', 'Control_Files', 
                                  'input_data', 'ex38c'])
    
    data_file = 'stratigraphy.csv'

    # Define keyword arguments of the system model
    num_years = 25
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    
    numberOfShaleLayers = 3
    
    # In the .csv file used, the center of the points is at x = 20000 km, y = 20000 km
    injX = 20000
    injY = 20000
    
    # Using x and y values outside of the range covered in the .csv file will 
    # cause LookupTableStratigraphy component to fail at the spatial interpolation.
    locX = 15000
    locY = 15000
    
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
    
    ares.add_par('logResPerm', value=-13.5, vary=False)
    ares.add_par('injRate', value=1, vary=False)
    # This value should not be too low relative to the locX and locY used
    ares.add_par('reservoirRadius', value=30000, vary=False)
    ares.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    
    # Link the unit thickness parameters to the LookupTableStratigraphy component
    for ob_nm in thickness_obs:
        ares.add_par_linked_to_obs(ob_nm, luts.linkobs[ob_nm])
    
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')
    
    msw = sm.add_component_model_object(MultisegmentedWellboreAI(
        name='msw', parent=sm))
    
    msw.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    
    # Link the unit thickness parameters to the LookupTableStratigraphy component
    for ob_nm in thickness_obs:
        msw.add_par_linked_to_obs(ob_nm, luts.linkobs[ob_nm])
    
    # Add parameters of Multisegmented Wellbore AI component
    msw.add_par('useDLmodel', value=1, vary=False)
    msw.add_par('wellRadius', value=0.035, vary=False)
    msw.add_par('logWellPerm', value=-12.5, vary=False)
    
    msw.add_obs('brine_aquifer1')
    msw.add_obs('brine_aquifer2')
    msw.add_obs('brine_atm')
    msw.add_obs('CO2_aquifer1')
    msw.add_obs('CO2_aquifer2')
    msw.add_obs('CO2_atm')
    msw.add_obs('mass_CO2_aquifer1')
    msw.add_obs('mass_CO2_aquifer2')
    
    msw.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    msw.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])
    
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
    
    print('Depths (m) to the bottom of each unit at x = {} km, y = {} km:'.format(
        locX / 1000, locY / 1000))
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
    
    CO2_aquifer1 = sm.collect_observations_as_time_series(msw, 'CO2_aquifer1')
    CO2_aquifer2 = sm.collect_observations_as_time_series(msw, 'CO2_aquifer2')
    CO2_atm = sm.collect_observations_as_time_series(msw, 'CO2_atm')
    
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
    
    mass_CO2_aquifer1 = sm.collect_observations_as_time_series(msw, 
                                                               'mass_CO2_aquifer1')
    mass_CO2_aquifer2 = sm.collect_observations_as_time_series(msw, 
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
    
