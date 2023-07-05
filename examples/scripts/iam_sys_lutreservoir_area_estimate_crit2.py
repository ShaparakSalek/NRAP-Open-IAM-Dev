# -*- coding: utf-8 -*-
"""
Example illustrates application of Area Estimate component to the reservoir data
provided through Lookup Table Reservoir component.

This example requires the additional FutureGen 2.0 data set.
FutureGen 2.0 data set can be downloaded from the following source:
https://edx.netl.doe.gov/dataset/phase-iii-nrap-open-iam

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/FutureGen2/1008_sims
"""
import os
import sys
import logging

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
import numpy as np
import matplotlib.pyplot as plt

from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, AreaEstimate)
from openiam.reservoir_data_interpolator import data_scatter_plot
from openiam.area_estimate_component import GRAVITY_CONSTANT


if __name__ == "__main__":

    __spec__ = None

    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'FutureGen2', '1008_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'fg1.csv'])):
        url = ''.join([
            'https://edx.netl.doe.gov/dataset/',
            'phase-iii-nrap-open-iam/resource/',
            '71aeb591-1609-430b-8392-4d75ee84750c \n'])
        msg = ''.join([
            '\nFutureGen 2.0 data set can be downloaded ',
            'from the following source:\n',
            url,
            'Check example description for more information.'])
        print(msg)

    # Define keyword arguments of the system model
    num_years = 1
    nt = 10
    time_array = 365.25*np.arange(0.0, num_years*nt+1)/nt
    sm_model_kwargs = {'time_array': time_array}  # time is given in days

    # Define stratigraphy
    num_aquifers = 4
    aquifer_thickness = [33.2, 84.1, 31.1, 61.6]
    num_shales = 5
    shale_thickness = [198.7, 74.4, 110.3, 118.9, 530.4]

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # Choose aquifer 2 as being aquifer of interest
    # Determine depth to the bottom of aquifer 2 using calculation from the reservoir depth
    # extracted from lookup tables
    lut_file_name = signature_data[1, -1]
    # column 2 has z coordinates data, column 4 has initial pressure
    file_data = np.loadtxt(os.path.join(file_directory, lut_file_name),
                      delimiter=',', skiprows=1, usecols=(2, 4))
    reservoir_depth = -file_data[:, 0]  # depth

    aquifer_depth = reservoir_depth - (shale_thickness[0] +
        aquifer_thickness[0] + shale_thickness[1])

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # The first row (except the last element) of the file contains names of the parameters
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    num_interpolators = 5

    # Create and add interpolators to the system model
    for ind in range(num_interpolators):
        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        interpolator1 = sm.add_interpolator(ReservoirDataInterpolator(
            name='int'+str(ind+1), parent=sm,
            header_file_dir=file_directory,
            time_file='time_points.csv',
            data_file='fg{}.csv'.format(ind+1),
            index=int(signature_data[ind+1, 0]),
            signature=signature), intr_family='reservoir')
        debug_msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(debug_msg)

    logging.debug('All interpolators are created')

    # Add reservoir component
    ltres = sm.add_component_model_object(
        LookupTableReservoir(name='ltres', parent=sm, intr_family='reservoir'))

    index_value = 4
    ltres.add_par('index', value=index_value, vary=False)

    output_directory = os.sep.join(['..', '..', 'output', 'area_estimate_test1'])
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Add gridded observations and observations to be linked
    for obs_nm in ['pressure', 'initial_pressure', 'CO2saturation']:
        # Add gridded observations of the reservoir component
        ltres.add_grid_obs(obs_nm, constr_type='array', output_dir=output_directory)
        # Add observations of reservoir component model to be used as input
        # of the AreaEstimate component
        ltres.add_obs_to_be_linked(obs_nm, obs_type='grid', constr_type='array')

    # Create Area Estimate component
    criteria = 2
    aest = sm.add_component_model_object(
        AreaEstimate(name='aest', parent=sm, criteria=criteria))

    # Add parameters of the component
    aquifer_brine_density = 1020.0
    # Reservoir is under pressured with respect to the aquifer
    res_brine_density = 1000.0

    aest.add_par('res_brine_density', value=res_brine_density, vary=False)
    aest.add_par('saturation_threshold', value=1.0e-10, vary=False)

    # Add additional keyword arguments of the component
    aquifer_pressure = 101325 + aquifer_depth*aquifer_brine_density*GRAVITY_CONSTANT
    aest.model_kwargs['aquifer_pressure'] = aquifer_pressure
    aest.model_kwargs['reservoir_depth'] = reservoir_depth
    # Aquifer depth is calculated based on the reservoir depth
    aest.model_kwargs['aquifer_depth'] = aquifer_depth

    # Add keyword arguments linked to gridded observations of the lookup table reservoir
    for kwarg_nm in ['pressure', 'initial_pressure', 'CO2saturation']:
        aest.add_kwarg_linked_to_obs(kwarg_nm, ltres.linkobs[kwarg_nm], obs_type='grid')

    # Add observations of the component
    for obs_nm in ['pressure_front', 'plume_extent', 'delineated_area']:
        aest.add_grid_obs(obs_nm, constr_type='array', output_dir=output_directory)

    # Run system model using current values of its parameters
    sm.forward()

    initial_pressure = sm.collect_gridded_observations_as_time_series(
        ltres, 'initial_pressure', output_directory,
        indices=[0], rlzn_number=0)[0] # only initial time point
    print('Initial pressure data shape:', initial_pressure.shape) # (3591,)
    pressure = sm.collect_gridded_observations_as_time_series(
        ltres, 'pressure', output_directory, rlzn_number=0) # all time points
    print('Pressure data shape:', pressure.shape)  # (11, 3591)
    saturation = sm.collect_gridded_observations_as_time_series(
        ltres, 'CO2saturation', output_directory, rlzn_number=0) # all time points
    print('CO2 saturation data shape:', saturation.shape)  # (11, 3591)

    pressure_front = sm.collect_gridded_observations_as_time_series(
        aest, 'pressure_front', output_directory, rlzn_number=0)
    plume_extent = sm.collect_gridded_observations_as_time_series(
        aest, 'plume_extent', output_directory, rlzn_number=0)
    delineated_area = sm.collect_gridded_observations_as_time_series(
        aest, 'delineated_area', output_directory, rlzn_number=0)

    # Close all plots if open
    plt.close('all')

    # Select time point to plot results
    t_ind = 5  # len(time_array)-1
    t_ind_p = time_array[t_ind]/365.25

    # Get x and y
    x = interpolator1.points[:, 0]
    y = interpolator1.points[:, 1]

    # Plot initial pressure
    fig1, ax1 = data_scatter_plot(
        x, y, initial_pressure, title='Initial pressure', cbar_label='Pa')

    # Plot pressure at the selected time point
    fig2, ax2 = data_scatter_plot(
        x, y, pressure[t_ind, :],
        title='Pressure at t = {} years'.format(t_ind_p),
        cbar_label='Pa')

    # Plot saturation at the selected time point
    xmin, xmax = 220000, 250000
    ymin, ymax = 4.390e+6, 4.430e+6
    pos_saturation_ind = np.where(saturation[t_ind, :]>0)[0]
    fig3, ax3 = data_scatter_plot(
        x[pos_saturation_ind], y[pos_saturation_ind],
        saturation[t_ind, :][pos_saturation_ind],
        title='CO2 saturation at t = {} years'.format(t_ind_p),
        cbar_label='Pa')
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(ymin, ymax)

    # Plot pressure change at the selected time point
    fig4, ax4 = data_scatter_plot(
        x, y, pressure[t_ind, :]-initial_pressure,
        title='Delta pressure at t = {} years'.format(t_ind_p),
        cbar_label='Pa')

    # Plot pressure front at the selected time point
    fig5, ax5 = data_scatter_plot(
        x, y, pressure_front[t_ind, :],
        title='Pressure front at t = {} years'.format(t_ind_p),
        cbar_label='', cmap='binary', vmin=0, vmax=1)
    ax5.set_xlim(xmin, xmax)
    ax5.set_ylim(ymin, ymax)

    # Plot plume front at the selected time point
    fig6, ax6 = data_scatter_plot(
        x, y, plume_extent[t_ind, :],
        title='Plume front at t = {} years'.format(t_ind_p),
        cbar_label='', cmap='binary', vmin=0, vmax=1)
    ax6.set_xlim(xmin, xmax)
    ax6.set_ylim(ymin, ymax)

    # Plot delineated area at the selected time point
    fig7, ax7 = data_scatter_plot(
        x, y, delineated_area[t_ind, :],
        title='Delineated area at t = {} years'.format(t_ind_p),
        cbar_label='', cmap='binary', vmin=0, vmax=1)
    ax7.set_xlim(xmin, xmax)
    ax7.set_ylim(ymin, ymax)
