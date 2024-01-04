"""
Example demonstrating the use of the LookupTableStratigraphy component. Here,
the component is used to interpolate unit thicknesses and depths at one
location.

Example of run:
$ python iam_sys_lutstrata.py
"""

import sys
import os
import logging
import numpy as np

from openiam.components.iam_base_classes import SystemModel
from openiam.components.stratigraphy_data_interpolator import StratigraphyDataInterpolator
from openiam.components.lookup_table_stratigraphy_component import LookupTableStratigraphy


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', 'Control_Files',
                                  'input_data', 'ex38a'])
    file_name = 'stratigraphy.csv'

    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    numberOfShaleLayers = 3

    # The location used has to fit within the domain shown in the file
    locX = 2000
    locY = 2000

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    intpr = sm.add_interpolator(StratigraphyDataInterpolator(
        name='int1', parent=sm,
        header_file_dir=file_directory,
        data_file=file_name),
        intr_family='stratigraphy')

    # Add stratigraphy component
    luts = sm.add_component_model_object(LookupTableStratigraphy(
        name='luts', parent=sm, intr_family='stratigraphy', locX=locX, locY=locY))

    luts.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)

    # Only use get_thickness_obs_names() and get_depth_obs_names() after the
    # numberOfShaleLayers parameter has been assigned.
    thickness_obs = luts.get_thickness_obs_names()

    # The thicness and depths observations are only produced in the first time
    # step, so use index=[0] when adding these observations.
    for ob_nm in thickness_obs:
        luts.add_obs(ob_nm, index=[0])

    depth_obs = luts.get_depth_obs_names()

    for ob_nm in depth_obs:
        luts.add_obs(ob_nm, index=[0])

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    # Collect thickness observations
    reservoirThickness = sm.collect_observations_as_time_series(
        luts, 'reservoirThickness', indices=[0])
    shale1Thickness = sm.collect_observations_as_time_series(
        luts, 'shale1Thickness', indices=[0])
    aquifer1Thickness = sm.collect_observations_as_time_series(
        luts, 'aquifer1Thickness', indices=[0])
    shale2Thickness = sm.collect_observations_as_time_series(
        luts, 'shale2Thickness', indices=[0])
    aquifer2Thickness = sm.collect_observations_as_time_series(
        luts, 'aquifer2Thickness', indices=[0])
    shale3Thickness = sm.collect_observations_as_time_series(
        luts, 'shale3Thickness', indices=[0])

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
    # named 'luts.reservoirDepth_#', where '_#' is the time step number. The
    # LookupTableStratigraphy component does not produce observations with
    # names that include the time step number.
    reservoirDepth = sm.collect_observations_as_time_series(
        luts, 'reservoirDepth', indices=[0])
    shale1Depth = sm.collect_observations_as_time_series(
        luts, 'shale1Depth', indices=[0])
    aquifer1Depth = sm.collect_observations_as_time_series(
        luts, 'aquifer1Depth', indices=[0])
    shale2Depth = sm.collect_observations_as_time_series(
        luts, 'shale2Depth', indices=[0])
    aquifer2Depth = sm.collect_observations_as_time_series(
        luts, 'aquifer2Depth', indices=[0])
    shale3Depth = sm.collect_observations_as_time_series(
        luts, 'shale3Depth', indices=[0])

    print('Depths (m) to the base of each unit:')
    print('    reservoirDepth: ', reservoirDepth)
    print('    shale1Depth: ', shale1Depth)
    print('    aquifer1Depth: ', aquifer1Depth)
    print('    shale2Depth: ', shale2Depth)
    print('    aquifer2Depth: ', aquifer2Depth)
    print('    shale3Depth: ', shale3Depth)

    reservoirMidDepth = sm.collect_observations_as_time_series(
        luts, 'reservoirMidDepth', indices=[0])
    shale1MidDepth = sm.collect_observations_as_time_series(
        luts, 'shale1MidDepth', indices=[0])
    aquifer1MidDepth = sm.collect_observations_as_time_series(
        luts, 'aquifer1MidDepth', indices=[0])
    shale2MidDepth = sm.collect_observations_as_time_series(
        luts, 'shale2MidDepth', indices=[0])
    aquifer2MidDepth = sm.collect_observations_as_time_series(
        luts, 'aquifer2MidDepth', indices=[0])
    shale3MidDepth = sm.collect_observations_as_time_series(
        luts, 'shale3MidDepth', indices=[0])

    print('Depths (m) to the middle of each unit:')
    print('    reservoirMidDepth: ', reservoirMidDepth)
    print('    shale1MidDepth: ', shale1MidDepth)
    print('    aquifer1MidDepth: ', aquifer1MidDepth)
    print('    shale2MidDepth: ', shale2MidDepth)
    print('    aquifer2MidDepth: ', aquifer2MidDepth)
    print('    shale3MidDepth: ', shale3MidDepth)

    reservoirTopDepth = sm.collect_observations_as_time_series(
        luts, 'reservoirTopDepth', indices=[0])
    shale1TopDepth = sm.collect_observations_as_time_series(
        luts, 'shale1TopDepth', indices=[0])
    aquifer1TopDepth = sm.collect_observations_as_time_series(
        luts, 'aquifer1TopDepth', indices=[0])
    shale2TopDepth = sm.collect_observations_as_time_series(
        luts, 'shale2TopDepth', indices=[0])
    aquifer2TopDepth = sm.collect_observations_as_time_series(
        luts, 'aquifer2TopDepth', indices=[0])
    shale3TopDepth = sm.collect_observations_as_time_series(
        luts, 'shale3TopDepth', indices=[0])

    print('Depths (m) to the top of each unit:')
    print('    reservoirTopDepth: ', reservoirTopDepth)
    print('    shale1TopDepth: ', shale1TopDepth)
    print('    aquifer1TopDepth: ', aquifer1TopDepth)
    print('    shale2TopDepth: ', shale2TopDepth)
    print('    aquifer2TopDepth: ', aquifer2TopDepth)
    print('    shale3TopDepth: ', shale3TopDepth)
