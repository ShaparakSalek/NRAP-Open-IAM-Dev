# -*- coding: utf-8 -*-
"""
Example illustrates application of Area Estimate component to the reservoir data
provided through Lookup Table Reservoir component.

This example requires the additional FutureGen 2.0 data set.
FutureGen 2.0 data set can be downloaded from the following source:
https://edx.netl.doe.gov/dataset/phase-iii-nrap-open-iam

The downloaded data set should be placed here:
    data/reservoir/lookuptables/FutureGen2/1008_sims
"""
import os
import sys
import logging

import numpy as np
import matplotlib.pyplot as plt

from openiam.components.iam_base_classes import SystemModel
from openiam.components.reservoir_data_interpolator import (
    data_scatter_plot, ReservoirDataInterpolator)
from openiam.components.lookup_table_reservoir_component import LookupTableReservoir
from openiam.components.area_estimate_component import AreaEstimate


if __name__ == "__main__":

    __spec__ = None

    file_directory = os.sep.join(['..', '..', 'data', 'reservoir',
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
    num_years = 20
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

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

    index_value = 5
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
    aest = sm.add_component_model_object(
        AreaEstimate(name='aest', parent=sm, criteria=1))

    # Add parameters of the component
    delta_pressure_threshold = 0.5e+6
    aest.add_par('delta_pressure_threshold', value=delta_pressure_threshold, vary=False)
    aest.add_par('saturation_threshold', value=1.0e-10, vary=False)

    # Add keyword arguments linked to gridded observations of the lookup table reservoir
    for kwarg_nm in ['pressure', 'initial_pressure', 'CO2saturation']:
        aest.add_kwarg_linked_to_obs(kwarg_nm, ltres.linkobs[kwarg_nm], obs_type='grid')

    # Add observations of the component
    for obs_nm in ['pressure_front', 'plume_extent', 'delineated_area']:
        aest.add_grid_obs(obs_nm, constr_type='array', output_dir=output_directory)
        aest.add_grid_obs('max_'+obs_nm, constr_type='array', output_dir=output_directory)

    # Run system model using current values of its parameters
    sm.forward()

    initial_pressure = sm.collect_gridded_observations_as_time_series(
        ltres, 'initial_pressure', output_directory,
        indices=[0], rlzn_number=0)[0] # only initial time point
    print('Initial pressure data shape:', initial_pressure.shape) # (3591,)
    pressure = sm.collect_gridded_observations_as_time_series(
        ltres, 'pressure', output_directory, rlzn_number=0) # all time points
    print('Pressure data shape:', pressure.shape)
    saturation = sm.collect_gridded_observations_as_time_series(
        ltres, 'CO2saturation', output_directory, rlzn_number=0) # all time points
    print('CO2 saturation data shape:', saturation.shape)

    pressure_front = sm.collect_gridded_observations_as_time_series(
        aest, 'pressure_front', output_directory, rlzn_number=0)
    plume_extent = sm.collect_gridded_observations_as_time_series(
        aest, 'plume_extent', output_directory, rlzn_number=0)
    delineated_area = sm.collect_gridded_observations_as_time_series(
        aest, 'delineated_area', output_directory, rlzn_number=0)
    max_delineated_area = sm.collect_gridded_observations_as_time_series(
        aest, 'max_delineated_area', output_directory, rlzn_number=0)

    # Get x and y
    plt.close('all')
    x = interpolator1.points[:, 0]
    y = interpolator1.points[:, 1]
    # Plot initial pressure
    fig1, ax1 = data_scatter_plot(
        x, y, initial_pressure, title='Initial pressure', cbar_label='Pa')

    # Final time point of simulation
    t_ind = len(time_array)-1

    # Plot pressure at the end of simulation
    fig2, ax2 = data_scatter_plot(
        x, y, pressure[t_ind, :],
        title='Pressure at t = {} years'.format(t_ind),
        cbar_label='Pa')

    # Plot saturation at the end of simulation
    xmin, xmax = 220000, 250000
    ymin, ymax = 4.390e+6, 4.430e+6
    pos_saturation_ind = np.where(saturation[t_ind, :]>0)[0]
    fig3, ax3 = data_scatter_plot(
        x[pos_saturation_ind], y[pos_saturation_ind],
        saturation[t_ind, :][pos_saturation_ind],
        title='CO2 saturation at t = {} years'.format(t_ind),
        cbar_label='Pa')
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(ymin, ymax)

    # Plot pressure change at the end of simulation
    fig4, ax4 = data_scatter_plot(
        x, y, pressure[t_ind, :]-initial_pressure,
        title='Delta pressure at t = {} years'.format(t_ind),
        cbar_label='Pa')

    # Plot pressure front at the end of simulation
    fig5, ax5 = data_scatter_plot(
        x, y, pressure_front[t_ind, :],
        title='Pressure front at t = {} years'.format(t_ind),
        cbar_label='', cmap='binary', vmin=0, vmax=1)
    ax5.set_xlim(xmin, xmax)
    ax5.set_ylim(ymin, ymax)

    # Plot plume front at the end of simulation
    fig6, ax6 = data_scatter_plot(
        x, y, plume_extent[t_ind, :],
        title='Plume front at t = {} years'.format(t_ind),
        cbar_label='', cmap='binary', vmin=0, vmax=1)
    ax6.set_xlim(xmin, xmax)
    ax6.set_ylim(ymin, ymax)

    # Plot delineated area at the end of simulation
    fig7a, ax7a = data_scatter_plot(
        x, y, delineated_area[t_ind, :],
        title='Delineated area at t = {} years'.format(t_ind),
        cbar_label='', cmap='binary', vmin=0, vmax=1)
    ax7a.set_xlim(xmin, xmax)
    ax7a.set_ylim(ymin, ymax)

    # Plot max delineated area at the end of simulation
    fig7b, ax7b = data_scatter_plot(
        x, y, max_delineated_area[t_ind, :],
        title='Maximum delineated area at t = {} years'.format(t_ind),
        cbar_label='', cmap='binary', vmin=0, vmax=1)
    ax7b.set_xlim(xmin, xmax)
    ax7b.set_ylim(ymin, ymax)

    # Plot pressure
    fig8, ax8 = plt.subplots(figsize=(8, 8))
    ax8.set_aspect('equal')
    ax8.set_title('Pressure')

    levels = np.linspace(0.8e+6, 2.0e+7, 30)
    im = ax8.tricontourf(x, y, pressure[t_ind, :], levels=levels)
    fig8.colorbar(im, label='Pa', ax=ax8)
    plt.show()

    # Combine two plots
    # Plot pressure change at the end of simulation
    fig9, ax9 = data_scatter_plot(
        x, y, pressure[t_ind, :]-initial_pressure,
        title='Delta pressure at t = {} years'.format(t_ind),
        cbar_label='Pa')
    levels = np.array([delta_pressure_threshold])
    ax9.tricontour(x, y, pressure[t_ind, :]-initial_pressure, levels=levels, colors='r')
