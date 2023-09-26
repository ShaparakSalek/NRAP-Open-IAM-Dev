# -*- coding: utf-8 -*-
"""
This script reads Decatur example data (available on EDX and tested in ORION)
in order to create lookup table in csv format acceptable in NRAP-Open-IAM
to create a possibility of linkage between NRAP-Open-IAM and ORION functionality,
as well as use the created data sets for separate analysis in NRAP-Open-IAM.

The script also reads and plots the seismic catalog data distributed with the
Decatur example in ORION. Both reservoir data and seismic catalog are needed
for the common type of analysis done in ORION.

ORION should be installed in environment through pip package
for all the imports to work properly or placed in the folder above root folder
for NRAP-Open-IAM.

The required data files decatur_pressure.hdf5 and decaturSeismic.hdf5 should be
placed in the folder source/components/reservoir/lookuptables/Decatur. The resulting
output lookup table and other files needed for the Lookup Table Reservoir component
will be placed there.

If you're an ORION user the data files can be downloaded from ORION workspace on EDX:
ask Chris Shermann (principal developer of ORION at shermann27@llnl.gov)
for details on how to access the examples and relevant data sets for ORION
If you're an NRAP-Open-IAM developer the original input and output data files
can be found on the developer's workspace on EDX "NRAP-Open-IAM Data Sharing":
ask the team lead for more information about this workspace and data.

The script was last tested on September 25th, 2023 with the most current version
of ORION (commit: f1c4bae38a4ca5c556848cff53c547bce06ac898) and NRAP-Open-IAM
(commit: 64e5802b4654896f350ec14c2af2ec0bf4d03c62).
"""

import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Check whether ORION is installed
try:
    from orion.utilities.hdf5_wrapper import hdf5_wrapper
    print(''.join(['Package orion is not installed. Trying to check if it is ',
                   'available in the folder above the root folder for NRAP-Open-IAM...']))
except ModuleNotFoundError:
    sys.path.insert(0, os.sep.join(['..', '..', '..', '..', '..', 'orion', 'src']))
    try:
        from orion.utilities.hdf5_wrapper import hdf5_wrapper
    except Exception as exc:
        raise ModuleNotFoundError('Package orion is not available on this system.') from exc

sys.path.insert(0, os.sep.join(['..', '..', '..', '..', 'source']))

from openiam import IAM_DIR


if __name__ == '__main__':
    # Setup input folder and required files names
    decatur_example_dir = os.sep.join([
        IAM_DIR, 'source', 'components', 'reservoir', 'lookuptables', 'Decatur'])
    figures_output_dir = os.sep.join([IAM_DIR, 'output', 'Decatur_example'])
    pressure_file_name = 'decatur_pressure.hdf5'
    seismic_file_name = 'decaturSeismic.hdf5'

    # Check that the data file needed for conversion is available
    if not os.path.exists(decatur_example_dir):
        raise FileNotFoundError('Folder {} is not found.'.format(decatur_example_dir))
    for file_nm in [pressure_file_name, seismic_file_name]:
        if not os.path.exists(os.sep.join([decatur_example_dir, file_nm])):
            raise FileNotFoundError(
                'Required file {} is not found in the folder {}.'.format(
                    file_nm, decatur_example_dir))

    # Create output folder for the figures
    if not os.path.exists(os.path.dirname(figures_output_dir)):
        os.mkdir(os.path.dirname(figures_output_dir))
    if not os.path.exists(figures_output_dir):
        os.mkdir(figures_output_dir)

    # Read pressure table file tested with ORION
    tmp = hdf5_wrapper(fname=os.sep.join([decatur_example_dir, pressure_file_name]))
    data = tmp.get_copy()

    for key in ['t', 'x', 'y', 'z', 'pressure', 'dpdt']:
        print('{} is in data = {}'.format(key, key in data))

    # Expected outcome
    # t is in data = True
    # x is in data = True
    # y is in data = True
    # z is in data = True
    # pressure is in data = True
    # dpdt is in data = True

    # Time points
    print('Number of time points', data['t'].shape)

    # x, y, and z coordinates
    nx = data['x'].shape[0]
    ny = data['y'].shape[0]
    nz = data['z'].shape[0]
    nt = data['t'].shape[0]

    x = data['x']
    y = data['y']
    z = data['z']

    # Transform to years from seconds expired from the epoch
    t = (data['t'] - min(data['t']))/(3600*24*365.25)
    dts = np.ediff1d(data['t'])

    init_time_point = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(data['t'][0]))
    print('Start of reservoir table model: ', init_time_point, data['t'][0])
    end_time_point = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(data['t'][-1]))
    print('End of reservoir table model: ', end_time_point, data['t'][-1])

    # Get representation of the reservoir data domain
    x_range = np.ptp(x)  # max-min
    y_range = np.ptp(y)
    z_range = np.ptp(z)
    x_min = np.min(x)
    y_min = np.min(y)
    z_min = np.min(z)
    x_max = np.max(x)
    y_max = np.max(y)
    z_max = np.max(z)
    dxs = np.ediff1d(x)
    dys = np.ediff1d(y)
    dzs = np.ediff1d(z)

    print('Domain size:')
    print('x in [{}, {}]'.format(x_min, x_max), 'Dx: {}'.format(x_range))
    print('y in [{}, {}]'.format(y_min, y_max), 'Dy: {}'.format(y_range))
    print('z in [{}, {}]'.format(z_min, z_max), 'Dz: {}'.format(z_range))

    print('Differences between consecutive points')
    print('dxs: ', dxs)
    print('dys: ', dys)
    print('dzs: ', dzs)

    xx1, yy1 = np.meshgrid(x, y, indexing='ij')
    xx2, zz2 = np.meshgrid(x, z, indexing='ij')
    print('Number of different x-coordinates', nx)
    print('Number of different y-coordinates', ny)
    print('Number of different z-coordinates', nz)

    # Check that the coordinates and temporal data were created correctly
    arrx, arry, arrz, arrt = np.meshgrid(x, y, z, t, indexing='ij')
    for ind1 in range(nx):
        for ind2 in range(ny):
            for ind3 in range(nz):
                for ind4 in range(nt):
                    if (arrx[ind1, ind2, ind3, ind4] != x[ind1]) or (
                                arry[ind1, ind2, ind3, ind4] != y[ind2]) or (
                                    arrz[ind1, ind2, ind3, ind4] != z[ind3]) or (
                                        arrt[ind1, ind2, ind3, ind4] != t[ind4]):
                        print('fail')

    # Setup pressure and saturation data for the lookup table
    pressure_data = np.zeros((nx*ny*nz, nt))
    for ind in range(nt):
        pressure_data[:, ind] = data['pressure'][:, :, :, ind].reshape((nx*ny*nz, ))

    saturation_data = np.zeros((nx*ny*nz, nt))
    x_data = arrx[:, :, :, 0].reshape((nx*ny*nz, 1))
    y_data = arry[:, :, :, 0].reshape((nx*ny*nz, 1))
    z_data = arrz[:, :, :, 0].reshape((nx*ny*nz, 1))

    all_data = np.concatenate((x_data, y_data, z_data, pressure_data, saturation_data), axis=1)
    header = ','.join(
        ['x', 'y', 'z'] + [
            'pressure_{}'.format(num) for num in range(1, nt+1)] + [
                'CO2saturation_{}'.format(num) for num in range(1, nt+1)])

    data_filename = os.sep.join([decatur_example_dir, "Reservoir_data_decatur.csv"])
    np.savetxt(data_filename, all_data, delimiter=",", header=header)

    # Save time points for csv lookup table
    t_file_name = os.sep.join([decatur_example_dir, 'time_points.csv'])
    with open(t_file_name, 'w') as f:
        f.write(','.join([str(val) for val in t])+'\n')

    signature_data = np.array([['index', 'filename'], [1, 'Reservoir_data_decatur.csv']])
    np.savetxt(os.sep.join([decatur_example_dir, 'parameters_and_filenames.csv']),
               signature_data, delimiter=",", fmt="%s")

    # Check that the data for the lookup table was created correctly
    for ind1 in range(nx):
        for ind2 in range(ny):
            for ind3 in range(nz):
                for ind4 in range(2):
                    if pressure_data[ind1*ny*nz+ind2*nz+ind3, ind4] != \
                            data['pressure'][ind1, ind2, ind3, ind4]:
                        print('fail')

    # Show range of x, y, z coordinates and domain of the example
    marker_size = 3
    to_show_domain = 1
    if to_show_domain:
        fig = plt.figure(1, figsize=(15, 6))
        ax1 = fig.add_subplot(121)
        ax1.plot(xx1.reshape((nx*ny, -1)), yy1.reshape((nx*ny, -1)), 'ko',
                 markersize=marker_size)
        ax1.set_aspect('equal')
        ax1.set_xlabel('x, [m]')
        ax1.set_ylabel('y, [m]')
        ax1.set_title('Top view')

        ax2 = fig.add_subplot(122)
        ax2.plot(xx2.reshape((nx*nz, -1)), zz2.reshape((nx*nz, -1)), 'ko',
                 markersize=marker_size)
        ax2.set_aspect('equal')
        ax2.set_xlabel('x, [m]')
        ax2.set_ylabel('depth, [m]')
        ax2.set_title('Side view')
        fig.suptitle('Decatur Example Domain')
        plt.savefig(os.sep.join([figures_output_dir, 'domain.png']),
                    dpi=300)

    to_process_seismic = 1
    if to_process_seismic:
        # Read seismic catalog
        s_temp = hdf5_wrapper(os.sep.join([decatur_example_dir, seismic_file_name]))
        s_data = s_temp.get_copy()

        print('Seismic data keys: ', s_data.keys())
        # Expected output
        # Seismic data keys:  dict_keys(['b_value', 'b_value_epoch', 'depth',
        # 'easting', 'epoch', 'event_id', 'latitude', 'longitude', 'magnitude',
        # 'northing', 'source_radius', 'stress_drop'])

        # time of event is in epoch (which is in seconds)
        eq_x = s_data['easting']
        eq_y = s_data['northing']
        magnitude = s_data['magnitude']
        latitude = s_data['latitude']
        longitude = s_data['longitude']

        to_show_eq_events = 1
        if to_show_eq_events:
            fig = plt.figure(2, figsize=(6, 8))
            ax1 = fig.add_subplot(111)
            domain_boundary = Rectangle((x_min, y_min), x_range, y_range,
                                        edgecolor='grey', facecolor='white')
            ax1.add_patch(domain_boundary)
            ax1.scatter(eq_x, eq_y, s=1+17*magnitude**2, c='k',
                        alpha=0.4)
            ax1.set_aspect('equal')
            ax1.set_xlabel('x, [m]')
            ax1.set_ylabel('y, [m]')
            fig.suptitle('Decatur seismic catalog')
            fig.tight_layout()
            plt.savefig(os.sep.join([figures_output_dir, 'seismic_catalog.png']),
                        dpi=300)
