"""
This example demonstrates the use of LookupTableStratigraphy component. This
component is used to produce unit thicknesses at multiple locations through
spatial interpolation of the data stored in a .csv file. At each location, the
unit thicknesses produced are linked to AnalyticalReservoir and OpenWellbore
components. For example, the AnalyticalReservoir component has parameters like
reservoirThickness, aquifer1Thickness, and shale3Thickness that are linked to
the observations of the LookupTableStratigraphy component. Additionally, the
wellTop parameter of the OpenWellbore component influences the critical pressure
calculated by the component, and the wellTop parameter is linked to the depth of
the aquifer receiving leakage.

Example of run:
$ python iam_sys_lutstrata_reservoir_openwell.py
"""

import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from openiam.components.iam_base_classes import SystemModel
from openiam.components.stratigraphy_data_interpolator import StratigraphyDataInterpolator
from openiam.components.lookup_table_stratigraphy_component import LookupTableStratigraphy
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.open_wellbore_component import OpenWellbore


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', 'Control_Files',
                                  'input_data', 'ex38a'])

    data_file = 'stratigraphy.csv'

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    numberOfShaleLayers = 3

    # In the .csv file used, the center of the points is at x = 3.5 km, y = 3.5 km
    injX = 0
    injY = 0

    # Using x and y values outside of the range covered in the .csv file will
    # cause LookupTableStratigraphy component to fail at the spatial interpolation.
    locX = [3000, 2000, 1000, 500, 250]
    locY = [0, 0, 0, 0, 0]
    loc_colors = ['mediumblue', 'yellowgreen', 'gold', 'orangered', 'darkred']

    # Plot the locations selected relative to the locations in the .csv file
    data = pd.read_csv(os.path.join(file_directory, data_file))

    all_x_vals = data.x
    all_y_vals = data.y

    plt.figure(1, figsize=(9, 7))

    ax = plt.gca()

    ax.plot(all_x_vals / 1000, all_y_vals / 1000, marker='o', color='k',
            markerfacecolor='none', markeredgewidth=2, markersize=4, linestyle='none',
            label='All Points From the File', zorder=1)

    for i, (loc_x, loc_y) in enumerate(zip(locX, locY)):
        ax.plot(loc_x / 1000, loc_y / 1000, marker='o', color=loc_colors[i],
                markerfacecolor='none', markeredgewidth=3, markersize=8, linestyle='none',
                label='Location {}'.format(i), zorder=3)

    ax.plot(injX / 1000, injY / 1000, marker='s', color='k', markerfacecolor='none',
            markeredgewidth=4, markersize=10, linestyle='none',
            label='Injection Site', zorder=2)

    ax.set_xlabel('Easting (km)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Northing (km)', fontsize=12, fontweight='bold')
    ax.set_title('Overview of the Domain', fontsize=12, fontweight='bold')
    plt.legend(fontsize=12)

    # The simulation will examine the leakage into this aquifer
    selected_aquifer_num = 2

    water_density_val = 1000
    brine_density_val = 1080
    grav_accel = 9.81

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    intpr = sm.add_interpolator(StratigraphyDataInterpolator(
        name='int1', parent=sm,
        header_file_dir=file_directory,
        data_file=data_file),
        intr_family='stratigraphy')

    luts = []
    ares = []
    ows = []
    for i, (loc_x, loc_y) in enumerate(zip(locX, locY)):
        # Add stratigraphy component
        luts.append(sm.add_component_model_object(LookupTableStratigraphy(
            name='luts{}'.format(i), parent=sm, intr_family='stratigraphy',
            locX=loc_x, locY=loc_y)))

        luts[i].add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)

        thickness_obs = luts[i].get_thickness_obs_names()

        for ob_nm in thickness_obs:
            luts[i].add_obs(ob_nm, index=[0])
            luts[i].add_obs_to_be_linked(ob_nm)

        depth_obs = luts[i].get_depth_obs_names()

        for ob_nm in depth_obs:
            luts[i].add_obs(ob_nm, index=[0])
            luts[i].add_obs_to_be_linked(ob_nm)

        ares.append(sm.add_component_model_object(AnalyticalReservoir(
            name='ares{}'.format(i), parent=sm, injX=injX, injY=injY,
            locX=loc_x, locY=loc_y)))

        ares[i].add_par('logResPerm', value=-12.5, vary=False)
        ares[i].add_par('injRate', value=0.05, vary=False)
        # This value should not be too low relative to the locX and locY used
        ares[i].add_par('reservoirRadius', value=10000, vary=False)
        ares[i].add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)

        # Link the unit thickness parameters to the LookupTableStratigraphy component
        for ob_nm in thickness_obs:
            ares[i].add_par_linked_to_obs(ob_nm, luts[i].linkobs[ob_nm])

        ares[i].add_obs('pressure')
        ares[i].add_obs('CO2saturation')
        ares[i].add_obs_to_be_linked('pressure')
        ares[i].add_obs_to_be_linked('CO2saturation')

        ows.append(sm.add_component_model_object(OpenWellbore(
            name='ow{}'.format(i), parent=sm, crit_pressure_approach=True,
            enforce_crit_pressure=False)))

        # Add parameters of multisegmented wellbore component
        ows[i].add_par('wellRadius', value=0.015, vary=False)
        ows[i].add_par('logReservoirTransmissivity', value=-10.1, vary=False)
        ows[i].add_par('logAquiferTransmissivity', value=-9.1, vary=False)
        ows[i].add_par('brineDensity', value=brine_density_val, vary=False)

        ows[i].add_par_linked_to_obs('reservoirDepth', luts[i].linkobs['shale1Depth'])
        ows[i].add_par_linked_to_obs('wellTop', luts[i].linkobs['aquifer{}Depth'.format(
            selected_aquifer_num)])

        ows[i].add_obs('brine_aquifer')
        ows[i].add_obs('CO2_aquifer')

        ows[i].add_kwarg_linked_to_obs('pressure', ares[i].linkobs['pressure'])
        ows[i].add_kwarg_linked_to_obs('CO2saturation', ares[i].linkobs['CO2saturation'])

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    # Plot the results
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': 10}
    plt.rc('font', **font)

    time_array_yrs = time_array / 365.25

    # This represents the type of unit extending from the depth in the same
    # index of unit_depths to the next index
    unit_types = ['Reservoir', 'Shale', 'Aquifer', 'Shale', 'Aquifer', 'Shale']
    unit_index = ['', ' 1', ' 1', ' 2', ' 2', ' 3']
    unit_colors = {'Reservoir': [0.33, 0.33, 0.33], 'Shale': 'olive', 'Aquifer': 'lightskyblue'}
    unit_x_range = {'Reservoir': [0, 8], 'Shale': [0, 7], 'Aquifer': [0, 8]}
    y_text_adjust = 10

    critP_list = []
    for i, (loc_x, loc_y) in enumerate(zip(locX, locY)):
        shale1Thickness = sm.collect_observations_as_time_series(luts[i], 'shale1Thickness',
                                                                 indices=[0])
        shale2Thickness = sm.collect_observations_as_time_series(luts[i], 'shale2Thickness',
                                                                 indices=[0])
        shale3Thickness = sm.collect_observations_as_time_series(luts[i], 'shale3Thickness',
                                                                 indices=[0])
        aquifer1Thickness = sm.collect_observations_as_time_series(luts[i], 'aquifer1Thickness',
                                                                   indices=[0])
        aquifer2Thickness = sm.collect_observations_as_time_series(luts[i], 'aquifer2Thickness',
                                                                   indices=[0])
        reservoirThickness = sm.collect_observations_as_time_series(luts[i], 'reservoirThickness',
                                                                    indices=[0])

        shale1Depth = sm.collect_observations_as_time_series(luts[i], 'shale1Depth',
                                                             indices=[0])
        shale2Depth = sm.collect_observations_as_time_series(luts[i], 'shale2Depth',
                                                             indices=[0])
        shale3Depth = sm.collect_observations_as_time_series(luts[i], 'shale3Depth',
                                                             indices=[0])
        aquifer1Depth = sm.collect_observations_as_time_series(luts[i], 'aquifer1Depth',
                                                               indices=[0])
        aquifer2Depth = sm.collect_observations_as_time_series(luts[i], 'aquifer2Depth',
                                                               indices=[0])

        if selected_aquifer_num == 1:
            aquifer_depth = aquifer1Depth[0]
        elif selected_aquifer_num == 2:
            aquifer_depth = aquifer2Depth[0]
        else:
            err_msg = ''.join(
                'The code needs to be edited to allow the calculated critical ',
                'pressure to reflect the selected aquifer number.')

            raise KeyError(err_msg)

        # This is how the OpenWellbore component calculates critical pressure
        critP = (aquifer_depth * grav_accel * water_density_val) + (
            brine_density_val * grav_accel * (shale1Depth[0] - aquifer_depth))

        critP_list.append(critP)

        # This list extends from the bottom of the reservoir to the surface
        unit_depths = [-(shale1Depth[0] + reservoirThickness[0]), -shale1Depth[0],
                       -aquifer1Depth[0], -shale2Depth[0], -aquifer2Depth[0],
                       -shale3Depth[0], 0]

        unit_thicknesses = [reservoirThickness[0], shale1Thickness[0],
                            aquifer1Thickness[0], shale2Thickness[0], aquifer2Thickness[0],
                            shale3Thickness[0]]

        plt.figure(2 + (i * 10), figsize=(4, 10))
        ax = plt.gca()
        for shaleRef in range(len(unit_depths) - 1):
            ax.plot(unit_x_range[unit_types[shaleRef]],
                    [unit_depths[shaleRef], unit_depths[shaleRef]],
                    color=unit_colors[unit_types[shaleRef]])

            ax.fill_between(
                unit_x_range[unit_types[shaleRef]],
                [unit_depths[shaleRef], unit_depths[shaleRef]],
                [unit_depths[shaleRef + 1], unit_depths[shaleRef + 1]],
                color=unit_colors[unit_types[shaleRef]], alpha=0.33)

            ax.text(unit_x_range[unit_types[shaleRef]][0], unit_depths[shaleRef] + y_text_adjust,
                    unit_types[shaleRef] + unit_index[shaleRef]
                    + ', Thickness: {:.2f} m'.format(unit_thicknesses[shaleRef]),
                    color=unit_colors[unit_types[shaleRef]], fontweight='bold')

        ax.set_xticks([])
        ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
        ax.set_title('Stratigraphy at Location {},\nx = {} km and y = {} km'.format(
            i, loc_x / 1000, loc_y / 1000), fontsize=12, fontweight='bold')

        # Plot pressure, CO2saturation, brine_aquifer, and CO2_aquifer
        pressure = sm.collect_observations_as_time_series(ares[i], 'pressure')
        CO2saturation = sm.collect_observations_as_time_series(ares[i], 'CO2saturation')

        brine_aquifer = sm.collect_observations_as_time_series(ows[i], 'brine_aquifer')
        CO2_aquifer = sm.collect_observations_as_time_series(ows[i], 'CO2_aquifer')

        plt.figure(3, figsize=(12, 8))

        ax = plt.subplot(2,1,1)
        ax.plot(time_array_yrs, pressure, color=loc_colors[i], linewidth=2,
                zorder=2)

        ax = plt.subplot(2,1,2)
        ax.plot(time_array_yrs, CO2saturation, color=loc_colors[i], linewidth=2,
                label='Location {}'.format(i))

        plt.figure(4, figsize=(12, 8))

        ax = plt.subplot(2,1,1)
        ax.plot(time_array_yrs, brine_aquifer, color=loc_colors[i], linewidth=2)

        ax = plt.subplot(2,1,2)
        ax.plot(time_array_yrs, CO2_aquifer, color=loc_colors[i], linewidth=2,
                label='Location {}'.format(i))


    plt.figure(3, figsize=(12, 8))

    ax = plt.subplot(2,1,1)

    # xlim = ax.get_xlim()
    for i, critP in enumerate(critP_list):
        ax.hlines(critP, min(time_array_yrs), max(time_array_yrs),
                  label='P$_{crit}$, ' + 'Location {}'.format(i),
                  color=loc_colors[i], linewidth=1, linestyle=':', zorder=1)
    # ax.set_xlim(xlim)

    ax.set_title('Reservoir Conditions', fontsize=12, fontweight='bold')
    ax.set_ylabel('Pressure (Pa)', fontsize=12, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
    plt.legend(fontsize=10)

    ax = plt.subplot(2,1,2)
    ax.set_ylabel('CO$_2$ Saturation', fontsize=12, fontweight='bold')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
    plt.legend(fontsize=10)

    plt.figure(4, figsize=(12, 8))

    ax = plt.subplot(2,1,1)
    ax.set_title('Leakage Rates to Aquifer {}'.format(selected_aquifer_num),
                  fontsize=12, fontweight='bold')
    ax.set_ylabel('Brine Leakage Rate (kg/s)',
                  fontsize=12, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')

    ax = plt.subplot(2,1,2)
    ax.set_ylabel('CO$_2$ Leakage Rate (kg/s)',
                  fontsize=12, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
    plt.legend(fontsize=10)
