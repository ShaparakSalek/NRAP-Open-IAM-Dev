# -*- coding: utf-8 -*-
"""
The stratigraphy data interpolator estimates unit thicknesses at (a) given
location(s) using user-prepared lookup tables. The calculations are based on
linear interpolation for irregular grids.
"""

import logging
import sys
import os
import h5py
import numpy as np
from numpy import matlib
from scipy.spatial import Delaunay
import pandas as pd

np.set_printoptions(threshold=sys.maxsize)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam.components.iam_gridded_observation import (interp_weights,
                                                            interpolate)
    from openiam.components.reservoir_data_interpolator import (check_file_format,
                                                                read_data_headers)
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


CANNOT_FIND_OBS_MSG = ''.join([
    'The LookupTableStratigraphy observation {} was not found within the observations ',
    'that are (1) constant over space and time or (2) those that are constant over ',
    'time but not space. LookupTableStratigraphy output cannot vary over time.',
    ' Check your input.'])


class StratigraphyDataInterpolator():
    """ NRAP-Open-IAM StratigraphyDataInterpolator class. """
    def __init__(self, name, parent, header_file_dir, data_file,
                 default_values=None, triangulation=None,
                 build_on_the_fly=False):
        """
        Constructor method of StratigraphyDataInterpolator

        :param name: name of interpolator
        :type name: str

        :param parent: name of the system model interpolator belongs to
        :type parent: SystemModel object

        :param header_file_dir: location (directory) of the stratigraphy
            data that will be used for interpolation
        :type header_file_dir: string

        :param data_file: name of *.csv or *.hdf5 file to read stratigraphy
            data from
        :type data_file: string

        :param default_values: dictionary of default values of observations
            which are not provided in the data_file. By default,
            it is assumed that data provided in the data_file is enough,
            thus, by default, the parameter is None.
            The values provided are assumed to be constant for all time steps
            and all spatial points.
        :type default_values: dict()

        :param triangulation: Delaunay tesselation in 2- or 3-d. By default,
            triangulation is None which means that for the given interpolator
            it has to be calculated based on the data available through
            the provided data_file.
        :type triangulation: scipy.spatial.Delaunay

        :param build_on_the_fly: flag variable indicating whether the data
            from the lookup table corresponding to the given interpolators
            will be read only if needed, e.g., at the first
            call of interpolator. By default, interpolator is created
            at the initialization of the corresponding StratigraphyDataInterpolator
            class object (value=False).
        :type build_on_the_fly: boolean

        :returns: StratigraphyDataInterpolator class object
        """
        # Setup attributes of the StratigraphyDataInterpolator class object
        self.name = name             # name of interpolator
        self.parent = parent         # a system model interpolator belongs to
        self.header_file_dir = header_file_dir  # path to the directory with simulation data files

        # data_file is a name of file with stratigraphy data
        self.hdf5_data_format, self.data_file = check_file_format(
            data_file, data_type='data')

        self.default_values = default_values  # values (if any) specified
        # for the interpolator observation

        # Setup default units of measurements for possible observations
        # (mainly for plotting purposes).
        # The dictionary can be updated later when we figure out whether any
        # additional observations are possible/needed; or dictionary can be updated
        # in script
        self.default_units = {'shaleThickness': 'm',
                              'aquiferThickness': 'm',
                              'reservoirThickness': 'm',
                              'depth': 'm'}

        self.output_bounds = {'shaleThickness': [1, 1600],
                              'aquiferThickness': [1, 1600],
                              'reservoirThickness': [1, 1600],
                              'depth': [-30000, -1]}

        # Save current triangulation
        self.triangulation = triangulation

        # Initialize a flag showing whether the data were read and interpolator was created
        self.created = False

        # Create models for each output type for each time step
        if not build_on_the_fly:
            self.create_data_interpolators(self.triangulation)

        # Log message about creating the interpolator
        msg = 'StratigraphyDataInterpolator created with name {name}'.format(name=name)
        logging.debug(msg)

    def process_csv_data_file(self):
        """
        Analyze the data provided in the csv data file and populate the corresponding
        attributes of the instance.
        """
        data = pd.read_csv(os.path.join(self.header_file_dir, self.data_file))

        # Initialize list of constant in space observations
        constant_obs = {}

        # Check whether z coordinate is present
        if self.data_headers[2] == 'z':
            all_names = self.data_headers[3:]
            delta_ind = 3
            constant_obs['depth'] = 2  # index of data in the csv file
        else:
            all_names = self.data_headers[2:]
            delta_ind = 2

        # Loop over all observation names in the header except the first two or three
        # (x, y and z if present)
        for ind, nm in enumerate(all_names):
            constant_obs[nm] = ind + delta_ind   # save index of the column

        for nm, ind in constant_obs.items():
            # -1 corresponds to the data constant in time but not in space
            self.data[nm] = {-1: data.iloc[:, ind].values}

        if self.default_values is not None:
            for nm, val in self.default_values.items():
                if nm in self.data:
                    msg = ''.join([
                        'Default value is specified for observation "{}" as ',
                        'argument of the StratigraphyDataInterpolator ',
                        'constructor method. Observation is also ',
                        'defined in the lookup data file "{}". Value provided ',
                        'to the constructor method will be used ',
                        'for calculations.']).format(nm, self.data_file)
                    logging.warning(msg)
                self.data[nm] = {-2: val}   # -2 corresponds to data constant both in space and time

        # Setup x, y points
        self.points = data.iloc[:, 0:2].values

        self.num_xy_points = len(self.points)

    def process_hdf5_data_file(self):
        """
        Analyze the data provided in the hdf5 data file and populate the corresponding
        attributes of the class instance.
        """
        # Initialize list of constant in space observations
        constant_obs = {}
        if 'z' in self.data_headers:
            constant_obs['depth'] = 'z'  # saving key where the data is kept

        coord_keys = ['x', 'y', 'z', 'ij', 'ijk']
        all_names = [nm for nm in self.data_headers if nm not in coord_keys]

        # Loop over all observation names in the header except coordinates
        # (x, y and z if present)
        for ind, nm in enumerate(all_names):
            constant_obs[nm] = nm   # save name of the data

        for nm, ind in constant_obs.items():
            # -1 corresponds to the data constant in time but not in space
            if nm != 'depth':
                self.data[nm] = {-1: nm}
            else:
                self.data[nm] = {-1: constant_obs[nm]}

        if self.default_values is not None:
            for nm, val in self.default_values.items():
                if nm in self.data:
                    msg = ''.join([
                        'Default value is specified for observation "{}" as ',
                        'argument of the StratigraphyDataInterpolator ',
                        'constructor method. Observation is also ',
                        'defined in the lookup data file "{}". Value provided ',
                        'to the constructor method will be used ',
                        'for calculations.']).format(nm, self.data_file)
                    logging.warning(msg)

                # -2 corresponds to data constant both in space and time
                self.data[nm] = {-2: val}

        data_file = os.path.join(self.header_file_dir, self.data_file)
        # Setup x, y (and z if needed) points
        with h5py.File(data_file, 'r') as hf:
            x = hf['x'][()]
            y = hf['y'][()]
            self.num_xy_points = len(x)

            self.points = np.zeros((self.num_xy_points, 2))

            self.points[:, 0] = x
            self.points[:, 1] = y

    def create_data_interpolators(self, triangulation=None):
        """
        Setup attributes storing applicable observation data.

        Create data attribute of the class object for use in the model method of the class.

        :param triangulation: Delaunay tesselation in 2-d. By default,
            triangulation is None which means that for the given interpolator
            it has to be calculated based on the data available through
            the provided data_file.
        :type triangulation: scipy.spatial.Delaunay
        """
        # Read data headers information
        self.data_headers = read_data_headers(self.hdf5_data_format,
                                              self.header_file_dir,
                                              self.data_file)

        # Setup data dictionary
        self.data = {}

        if not self.hdf5_data_format:     # csv data file format
            self.process_csv_data_file()
        else:                             # hdf5 data file format
            self.process_hdf5_data_file()

        # Determine triangulation for the data points
        if triangulation is None:
            self.triangulation = Delaunay(self.points)

        # Record that the data is read and the interpolator is created
        self.created = True

    def data_units(self, obs_nm):
        """
        Return measurement units of the provided observations.

        :param obs_nm: name of the observation for which units of measurements
            are to be provided.
        :type obs_nm: str

        :returns: units of measurements as str. For observation not
            in the list of possible known names, the string 'Unknown' is returned.
        """
        if 'shale' in obs_nm and 'Thickness' in obs_nm:
            obs_type = 'shaleThickness'
        elif 'aquifer' in obs_nm and 'Thickness' in obs_nm:
            obs_type = 'aquiferThickness'
        else:
            obs_type = obs_nm

        if obs_type in self.default_units:
            return self.default_units[obs_type]

        return 'Unknown'

    def check_output_bounds(self, obs_nm, val):
        """
        This function checks if an observation value falls outside the boundaries
        of the observation type. If it does, the observation is set to the
        closest boundary.
        """
        if 'shale' in obs_nm and 'Thickness' in obs_nm:
            obs_type = 'shaleThickness'
        elif 'aquifer' in obs_nm and 'Thickness' in obs_nm:
            obs_type = 'aquiferThickness'
        else:
            obs_type = obs_nm

        unit = self.default_units[obs_type]

        warning_msg_pt1 = ''.join([
            'When using spatial interpolation on the stratigraphy data provided, ',
            'the observation {} was found to have a value of {} {}. '])

        warning_msg_pt2 = ''.join([
            'This value falls outside of the boundaries allowed for that ',
            'observation type, {} to {} {}. ']).format(
                self.output_bounds[obs_type][0], self.output_bounds[obs_type][1],
                unit)

        warning_msg_pt3 = 'The value will be set to {} {}. Check your input.'

        try:
            if val < self.output_bounds[obs_type][0]:
                warning_msg = warning_msg_pt1.format(obs_nm, val, unit) + \
                    warning_msg_pt2 + warning_msg_pt3.format(
                        self.output_bounds[obs_type][0], unit)
                logging.warning(warning_msg)

                val = np.array([self.output_bounds[obs_type][0]])

            elif val > self.output_bounds[obs_type][1]:
                warning_msg = warning_msg_pt1.format(obs_nm, val, unit) + \
                    warning_msg_pt2 + warning_msg_pt3.format(
                        self.output_bounds[obs_type][1], unit)

                logging.warning(warning_msg)

                val = np.array([self.output_bounds[obs_type][1]])

        except:
            for ind, value in enumerate(val):
                if value < self.output_bounds[obs_type][0]:
                    warning_msg = warning_msg_pt1.format(obs_nm, value, unit) + \
                        warning_msg_pt2 + warning_msg_pt3.format(
                            self.output_bounds[obs_type][0], unit)
                    logging.warning(warning_msg)

                    val[ind] = np.array([self.output_bounds[obs_type][0]])

                elif value > self.output_bounds[obs_type][1]:
                    warning_msg = warning_msg_pt1.format(obs_nm, value, unit) + \
                        warning_msg_pt2 + warning_msg_pt3.format(
                            self.output_bounds[obs_type][1], unit)

                    logging.warning(warning_msg)

                    val[ind] = np.array([self.output_bounds[obs_type][1]])

        return val

    def calculate_csv_based_output(self, vtx=None, wts=None):
        """
        Calculate outputs using approach based on csv files.
        """
        # Create dictionary of output
        out = dict()

        # Check whether a single point is requested
        if vtx is not None and wts is not None:
            # Check whether weights are reasonable; if they are not, then it means
            # that the point at which we need to interpolate is outside the data range
            if np.any(abs(wts) > 1.0+1.0e-8):
                msg = ''.join(['Input array of weights {} has invalid ',
                               'values: some or all values are ',
                               'greater than 1.']).format(wts)
                logging.warning(msg)

            # Cycle over all observations
            for nm, data_item in self.data.items():
                if -2 in data_item:     # obs is constant in space and time
                    val = np.array([data_item[-2]])
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val

                elif -1 in data_item:   # obs is constant in time but not space
                    val = interpolate(data_item[-1], vtx, wts)
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val

                else:                   # obs varies in space and time
                    logging.warning(CANNOT_FIND_OBS_MSG.format(nm))

        else:  # return data corresponding to all points in the data file
            for nm, data_item in self.data.items():
                if -2 in data_item:     # obs is constant in space and time
                    val = data_item[-2]
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val * np.ones(self.num_xy_points)

                elif -1 in data_item:   # obs is constant in time but not space
                    val = data_item[-1]
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val

                else:                   # obs varies in space and time
                    logging.warning(CANNOT_FIND_OBS_MSG.format(nm))

        return out

    def calculate_hdf5_based_output(self, vtx=None, wts=None):
        """
        Calculate outputs using approach based on csv files.
        """
        # Create dictionary of output
        out = dict()

        # Get path to data file
        data_file = os.path.join(self.header_file_dir, self.data_file)

        # Check whether a single point is requested
        if vtx is not None and wts is not None:
            # Check whether weights are reasonable; if they are not, then it means
            # that the point at which we need to interpolate is outside the data range
            if np.any(abs(wts) > 1.0+1.0e-8):
                msg = ''.join(['Input array of weights {} has invalid ',
                               'values: some or all values are ',
                               'greater than 1.']).format(wts)
                logging.warning(msg)

            # Cycle over all observations
            for nm, data_item in self.data.items():
                if -2 in data_item:     # obs is constant in space and time
                    val = np.array([data_item[-2]])
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val

                elif -1 in data_item:   # obs is constant in time but not space
                    with h5py.File(data_file, 'r') as hf:
                        data = hf[data_item[-1]][()]
                    val = interpolate(data, vtx, wts)
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val

                else:                   # obs varies in space and time
                    logging.warning(CANNOT_FIND_OBS_MSG.format(nm))

        else:  # return data corresponding to all points in the data file
            for nm, data_item in self.data.items():
                if -2 in data_item:     # obs is constant in space and time
                    val = data_item[-2]
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val * np.ones(self.num_xy_points)

                elif -1 in data_item:   # obs is constant in time but not space
                    with h5py.File(data_file, 'r') as hf:
                        data = hf[data_item[-1]][()]
                    val = data
                    val = self.check_output_bounds(nm, val)
                    out[nm] = val

                else:                   # obs varies in space and time
                    logging.warning(CANNOT_FIND_OBS_MSG.format(nm))

        return out

    def __call__(self, vtx=None, wts=None):
        """
        Return observation data at the point of interest.

        :param vtx: array of indices of simplex vertices which enclose
            the point of interest; array should have a shape (1,3); indices
            do not exceed the number of data points
        :type vtx: numpy.array of int

        :param wts: array of weights of data at simplex vertices which enclose
            the point of interest; array should have a shape (1,3); weights
            values should be between 0 and 1, and sum up to 1
        :type wts: numpy.array of floats

        :returns: out - dictionary of observations; keys depend on the data
            provided in the data files. Possible keys are (here, N is the number
            of shale layers):
            ['shale1Thicness', 'shale2Thickness', ..., 'shaleNThickness',
             'aquifer1Thickness', 'aquifer2Thickness', ..., 'aquiferN-1Thickness',
             'reservoirThickness'].
        """
        if not self.created:
            self.create_data_interpolators(self.triangulation)

        if self.hdf5_data_format:
            out = self.calculate_hdf5_based_output(vtx=vtx, wts=wts)
        else:
            out = self.calculate_csv_based_output(vtx=vtx, wts=wts)

        return out


def test_case1():
    try:
        from openiam.components.iam_base_classes import SystemModel
    except ImportError as err:
        print('Unable to load NRAP-Open-IAM base classes module: {}'.format(err))

    logging.basicConfig(level=logging.WARNING)

    min_val = 0
    max_val = 20000
    interval = 1000

    coord_vals =  np.arange(min_val, max_val + interval, interval)

    x_vals = matlib.repmat(coord_vals, 1, len(coord_vals))
    x_vals = x_vals.reshape(-1, 1)[:, 0]

    y_vals = []
    y_val = 0
    for i in range(len(coord_vals)):
        for j in range(len(coord_vals)):
             y_vals.append(y_val)

        y_val += interval

    y_vals = np.array(y_vals)

    reservoirThickness = np.zeros(len(x_vals))
    shale1Thickness = np.zeros(len(x_vals))
    aquifer1Thickness = np.zeros(len(x_vals))
    shale2Thickness = np.zeros(len(x_vals))
    aquifer2Thickness = np.zeros(len(x_vals))
    shale3Thickness = np.zeros(len(x_vals))
    for i, (x_val, y_val) in enumerate(zip(x_vals, y_vals)):
        reservoirThickness[i] = 40 + (0.001 * x_val) + (0.0005 * y_val)
        shale1Thickness[i] = 170 + (0.005 * x_val) + (0.00075 * y_val)
        aquifer1Thickness[i] = 35 + (0.0025 * x_val) + (0.00075 * y_val)
        shale2Thickness[i] = 120 + (0.00125 * x_val) + (0.001 * y_val)
        aquifer2Thickness[i] = 75 - (0.0005 * x_val) - (0.00125 * y_val)
        shale3Thickness[i] = 220 - (0.00125 * x_val) - (0.005 * y_val)

    # In case the dips above are changed, prevent the units from becoming too thin
    reservoirThickness[reservoirThickness < 10] = 10
    shale1Thickness[shale1Thickness < 25] = 25
    aquifer1Thickness[aquifer1Thickness < 10] = 10
    shale2Thickness[shale2Thickness < 25] = 25
    aquifer2Thickness[aquifer2Thickness < 10] = 10
    shale3Thickness[shale3Thickness < 25] = 25

    # A .csv file will be written to this directory
    file_directory = os.sep.join([os.getcwd(), '..', '..', 'examples', 'scripts',
                                  'stratigraphy_data_interpolator'])

    if not os.path.exists(file_directory):
        os.mkdir(file_directory)

    file_name = 'stratigraphy.csv'
    file_path = os.path.join(file_directory, file_name)

    data = {'x': x_vals,
            'y': y_vals,
            'reservoirThickness': reservoirThickness,
            'shale1Thickness': shale1Thickness,
            'aquifer1Thickness': aquifer1Thickness,
            'shale2Thickness': shale2Thickness,
            'aquifer2Thickness': aquifer2Thickness,
            'shale3Thickness': shale3Thickness}

    strat_data = pd.DataFrame(data=data)

    strat_data.to_csv(file_path, index=False)

    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # The location used must fit within the min_val and max_val defined above
    locX = 12000
    locY = 10000

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add stratigraphy data interpolator
    intpr = sm.add_interpolator(StratigraphyDataInterpolator(
        name='int1', parent=sm,
        header_file_dir=file_directory,
        data_file=file_name),
        intr_family='stratigraphy')

    # Setup location of interest, must fit within the x and y values defined above
    locX, locY = [12000.0, 10000.0]

    # Calculate weights of the location of interest
    vertices, weights = interp_weights(intpr.triangulation, np.array([[locX, locY]]))

    keys = ['reservoirThickness', 'shale1Thickness', 'aquifer1Thickness',
            'shale2Thickness', 'aquifer2Thickness', 'shale3Thickness']

    out = intpr(vertices, weights)

    print('At x = {} m and y = {} m, the unit thicknesses are: '.format(locX, locY))

    for key in keys:
        print('    ' + key + ': ', out[key][0], ' m')


def test_case2():
    pass


if __name__ == "__main__":

    test_case1()
