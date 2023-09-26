"""
Methods included allow setting up lookup table reservoir component and provide
as ouput the pressure data in h5 file format accepted by ORION.

The included method test_Decatur_scenario uses Decatur example lookup table
data for the testing purposes. The three required data files
(parameters_and_filenames.csv, Reservoir_data_decatur.csv, and time_points.csv)
should be placed in the folder source/components/reservoir/lookuptables/Decatur.

If you're an NRAP-Open-IAM developer the needed data files
can be found on the developer's workspace on EDX "NRAP-Open-IAM Data Sharing":
ask the team lead for more information about this workspace and data.
"""

import sys
import os
import logging
import time
import numpy as np
import matplotlib.pyplot as plt
from orion.utilities import hdf5_wrapper
sys.path.insert(0, os.sep.join(['..', '..', '..', 'source']))

from openiam import IAM_DIR, SystemModel, ReservoirDataInterpolator, LookupTableReservoir
from matk import pyDOE

DEFAULT_OUTPUT_DIR = os.sep.join([IAM_DIR, 'output', 'orion_data'])

INITIAL_DATE_TEMPLATE = {key: 0 for key in ['hour', 'minute', 'second']}


def create_lutr_comp_setup(time_points, coord_x, coord_y, coord_z,
                           reservoir_pars, initial_true_date,
                           file_directory=None,
                           time_file='time_points.csv',
                           parameter_filename='parameters_and_filenames.csv',
                           output_directory=DEFAULT_OUTPUT_DIR,
                           orion_data_file='orion_data.hdf5'):
    """
    Parameters
    ----------
    time_points : array-like, e.g. list or numpy.array
        array of "relative" time points at which pressure and saturation data is
        requested in seconds; time_points[0] = 0
    coord_x : array-like, e.g. list or numpy.array
        unique x-coordinates of points at which reservoir data is to be produced
    coord_y : array-like, e.g. list or numpy.array
        unique y-coordinates of points at which reservoir data is to be produced
    coord_z : array-like, e.g. list or numpy.array
        unique z-coordinates of points at which reservoir data is to be produced
    reservoir_pars : dict
        Parameters names and values that are associated with a used lookup table file
    initial_true_date : true (real world) date corresponding to t=0 defined as dictionary
        with keys 'year', 'month', 'day', 'hour', 'minute', 'second'.
        Keys 'hour', 'minute', and 'second' are not required. If not provided
        they are assumed to be zero.
    file_directory : str, optional
        Path to the folder containing lookup table data file. The default is None.
    time_file : str, optional
        Name of the file containing time points associated with the lookup table
        data file. The default is 'time_points.csv'.
    parameter_filename : TYPE, optional
        Name of the file containing information about parameters and files names
        associated with the lookup table data set. The default is 'parameters_and_filenames.csv'.
    output_directory : str, optional
        Path to the folder where the produced pressure tables for ORION will
        be saved. The default is DEFAULT_OUTPUT_DIR.
    orion_data_file : str, optional
        Name of the file that will keep the produced pressure table for ORION.
        The default is 'orion_data.hdf5'.

    Returns
    -------
    None.

    """
    # Define keyword arguments of the system model
    sm_model_kwargs = {'time_array': time_points/86400.0}   # time is provided in days to the system model
    nt = len(time_points)

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, parameter_filename),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header

    # Create and add interpolators to the system model
    for ind in range(num_interpolators):
        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        interpltr = sm.add_interpolator(ReservoirDataInterpolator(
            name='int'+str(ind+1), parent=sm,
            header_file_dir=file_directory,
            time_file=time_file,
            data_file=signature_data[ind+1, -1],
            index=int(signature_data[ind+1, 0]),
            signature=signature,
            interp_2d=False), intr_family='reservoir')
        debug_msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(debug_msg)

    logging.debug('All interpolators are created')

    # Determine if reservoir parameters are defined by index
    index = reservoir_pars.pop('index', -1)

    # We assume we deal with 3d reservoir data
    nx = len(coord_x)
    ny = len(coord_y)
    nz = len(coord_z)
    num_points = nx*ny*nz
    print('Total number of requested points:', num_points)

    x_grid, y_grid, z_grid = np.meshgrid(coord_x, coord_y, coord_z, indexing='ij')

    x_grid = x_grid.reshape((num_points, ))
    y_grid = y_grid.reshape((num_points, ))
    z_grid = z_grid.reshape((num_points, ))

    ltres = sm.add_component_model_object(
        LookupTableReservoir(name='ltres', parent=sm,
                            locX=x_grid, locY=y_grid, locZ=z_grid,
                            intr_family='reservoir',
                            interp_2d=False))

    if index == -1:
        for par_name, par_val in reservoir_pars.items():
            ltres.add_par(par_name, value=par_val, vary=False)
    else:
        ltres.add_par('index', value=index, vary=False)

    for obs_name in ['pressure', 'CO2saturation']:
        ltres.add_grid_obs(obs_name, constr_type='array', output_dir=output_directory)

    # Run simulation
    sm.forward()

    # Collect observations
    output = {}
    for obs_name in ['pressure', 'CO2saturation']:
        # output[obs_name] will be of size num_time_points x num_points
        # So we'll need to transpose and reshape it before saving
        output[obs_name] = sm.collect_gridded_observations_as_time_series(
            ltres, obs_name, output_directory,
            indices=None, rlzn_number=0)
        output[obs_name] = output[obs_name].T.reshape((nx, ny, nz, nt))

    # Convert initial time point into epoch seconds
    initial_date = INITIAL_DATE_TEMPLATE.copy()
    initial_date.update(initial_true_date)
    t0_sec = time.mktime((initial_date['year'],
                          initial_date['month'],
                          initial_date['day'],
                          initial_date['hour'],
                          initial_date['minute'],
                          initial_date['second'], 0, 0, 0))
    print(t0_sec)
    save_t = np.array(time_points) + t0_sec
    print('First time point: ', save_t[0])
    print('Last time point: ', save_t[-1])

    # Check dimensions
    print('x shape: ', coord_x.shape)
    print('y shape: ', coord_y.shape)
    print('z shape: ', coord_z.shape)
    print('t shape: ', save_t.shape)
    print('pressure shape: ', output['pressure'].shape)

    if not os.path.exists(os.path.dirname(output_directory)):
        os.mkdir(os.path.dirname(output_directory))
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Save file for ORION
    orion_output_file = os.sep.join([output_directory, orion_data_file])
    with hdf5_wrapper.hdf5_wrapper(orion_output_file, mode='w') as data:
        data['x'] = coord_x
        data['y'] = coord_y
        data['z'] = coord_z
        data['t'] = save_t
        data['p'] = output['pressure']


def test_Decatur_scenario(test=2):
    file_directory = os.sep.join(['..', '..', '..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Decatur'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_decatur.csv'])):
        msg = ''.join(['\nRequired data set is not found.',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Initial time point in pressure table model for ORION Decatur example
    t0 = 1321516800.0
    # Last time point in pressure table model for ORION Decatur example
    tn = 1419732773.4375
    # Difference between tn and t0 is 98215973.4375 seconds. For approximation
    # we take 9.8e+7
    T = 98000000.0  # seconds converted to years gives 3.1; to days gives 1134.25
    dt = 1.0e+6  # seconds; difference between consecutive time points for simulation

    # time is provided in seconds to the method and converted to days
    # for input to the system model
    time_points = np.linspace(0, T, num=99)

    if test == 1:
        # Original grid
        coord_x = np.array([335580.5, 335622.8968254, 335673.29365079, 335723.69047619,
                            335774.08730159, 335824.48412698, 335874.88095238, 335925.27777778,
                            335975.67460317, 336026.07142857, 336076.46825397, 336126.86507937,
                            336177.26190476, 336227.65873016, 336278.05555556, 336328.45238095,
                            336378.84920635, 336429.24603175, 336479.64285714, 336530.03968254,
                            336580.43650794, 336630.83333333, 336681.23015873, 336731.62698413,
                            336782.02380952, 336832.42063492, 336882.81746032, 336933.21428571,
                            336983.61111111, 337034.00793651, 337084.4047619, 337134.8015873 ,
                            337185.1984127, 337235.5952381, 337285.99206349, 337336.38888889,
                            337386.78571429, 337437.18253968, 337487.57936508, 337537.97619048,
                            337588.37301587, 337638.76984127, 337689.16666667, 337739.56349206,
                            337789.96031746, 337840.35714286, 337890.75396825, 337941.15079365,
                            337991.54761905, 338041.94444444, 338092.34126984, 338142.73809524,
                            338193.13492063, 338243.53174603, 338293.92857143, 338344.32539683,
                            338394.72222222, 338445.11904762, 338495.51587302, 338545.91269841,
                            338596.30952381, 338646.70634921, 338697.5, 338747.5])[0:-1:3]
        coord_y = np.array([4415104.5, 4415054.95454545, 4415105.40909091, 4415155.86363636,
                            4415206.31818182, 4415256.77272727, 4415307.22727273, 4415357.68181818,
                            4415408.13636364, 4415458.59090909, 4415509.04545455, 4415559.5,
                            4415609.95454545, 4415660.40909091, 4415710.86363636, 4415761.31818182,
                            4415811.77272727, 4415862.22727273, 4415912.68181818, 4415963.13636364,
                            4416013.59090909, 4416064.04545455, 4416114.5, 4416164.95454545,
                            4416215.40909091, 4416265.86363636, 4416316.31818182, 4416366.77272727,
                            4416417.22727273, 4416467.68181818, 4416518.13636364, 4416568.59090909,
                            4416619.04545455, 4416669.5, 4416719.95454545, 4416770.40909091,
                            4416820.86363636, 4416871.31818182, 4416921.77272727, 4416972.22727273,
                            4417022.68181818, 4417073.13636364, 4417123.59090909, 4417174.04545455,
                            4417224.5, 4417274.95454545, 4417325.40909091, 4417375.86363636,
                            4417426.31818182, 4417476.77272727, 4417527.22727273, 4417577.68181818,
                            4417628.13636364, 4417678.59090909, 4417729.5, 4417779.5])[0:-1:3]
        coord_z = np.array([-2968.75, -2917.515625, -2866.28125, -2815.046875, -2763.8125,
                            -2712.578125, -2661.34375, -2610.109375, -2558.875, -2507.640625,
                            -2456.40625, -2405.171875, -2353.9375, -2302.703125, -2251.46875,
                            -2200.234375, -2149.])
        reservoir_pars = {'index': 1}

        create_lutr_comp_setup(
            time_points, coord_x, coord_y, coord_z, reservoir_pars,
            initial_true_date={'year': 2011, 'month': 11, 'day': 17},  # November 17th, 2011
            file_directory=file_directory,
            time_file='time_points.csv',
            parameter_filename='parameters_and_filenames.csv',
            output_directory=os.sep.join([IAM_DIR, 'output', 'Decatur_example']),
            orion_data_file = 'orion_data.hdf5')
    else:
        # Reduced and less dense grid
        coord_x = np.linspace(335600.0, 338700.0, num=32)
        coord_y = np.linspace(4415100.0, 4417700.0, num=27)
        coord_z = np.array([-2968.75, -2917.515625, -2866.28125, -2815.046875, -2763.8125,
                            -2712.578125, -2661.34375, -2610.109375, -2558.875, -2507.640625,
                            -2456.40625, -2405.171875, -2353.9375, -2302.703125, -2251.46875,
                            -2200.234375, -2149.])

        reservoir_pars = {'index': 1}

        create_lutr_comp_setup(
            time_points, coord_x, coord_y, coord_z, reservoir_pars,
            initial_true_date={'year': 2011, 'month': 11, 'day': 17},  # November 17th, 2011
            file_directory=file_directory,
            time_file='time_points.csv',
            parameter_filename='parameters_and_filenames.csv',
            output_directory=os.sep.join([IAM_DIR, 'output', 'Decatur_example']),
            orion_data_file = 'orion_data_selected.hdf5')


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)

    test = 2

    test_Decatur_scenario(test=test)
