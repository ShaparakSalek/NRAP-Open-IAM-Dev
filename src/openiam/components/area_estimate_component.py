# -*- coding: utf-8 -*-
import logging
import sys
import os
import numpy as np

try:
    import openiam.components.iam_base_classes as iam_bc
except ImportError as err:
    print('Unable to load NRAP-Open-IAM base classes module: {}'.format(err))


AE_OBSERVATIONS = ['pressure_front', 'plume_extent', 'delineated_area']

PARAMETERS_GROUPS =  {
    1: ['delta_pressure_threshold', 'saturation_threshold'],
    2: ['saturation_threshold', 'aquifer_pressure',
        'res_brine_density', 'reservoir_depth', 'aquifer_depth'],
    3: ['saturation_threshold', 'res_brine_density', 'aq_brine_density',
        'reservoir_depth', 'aquifer_depth']}

GRAVITY_CONSTANT = 9.8


def simple_delineation(p, pressure, saturation, initial_pressure):
    """
    Calculate pressure front, CO2 plume extent and delineated area.

    :param p: input parameters
    :type p: dict()

    :param pressure: pressure in the reservoir
    :type pressure: numpy.array of shape (N, ) where N is a number of points

    :param saturation: CO2 saturation in the reservoir
    :type saturation: numpy.array of shape (N, ) where N is a number of points

    :param initial_pressure: pressure in the reservoir before injection
    :type initial_pressure: numpy.array of shape (N, ) where N is a number of points

    :returns:
        out: Dictionary with keys 'pressure_front', 'plume_extent', and 'delineated_area'.

    """
    # Calculate difference in the initial and current pressures
    delta_pressure = pressure - initial_pressure

    # Calculate pressure front
    pressure_front = np.zeros(pressure.shape)
    pressure_front[np.where(delta_pressure > p['delta_pressure_threshold'])] = 1

    # Calculate CO2 plume extent
    plume_extent = np.zeros(saturation.shape)
    plume_extent[np.where(saturation > p['saturation_threshold'])] = 1

    delineated_area = np.logical_or(pressure_front, plume_extent)

    # Setup dictionary of outputs
    out = {'pressure_front': pressure_front,
           'plume_extent': plume_extent,
           'delineated_area': delineated_area}

    return out

def under_pressurized_case_delineation(p, pressure, saturation, initial_pressure):
    """
    Calculate pressure front, CO2 plume extent and delineated area

    :param p: input parameters
    :type p: dict()

    :param pressure: pressure in the reservoir
    :type pressure: numpy.array of shape (N, ) where N is a number of points

    :param saturation: CO2 saturation in the reservoir
    :type saturation: numpy.array of shape (N, ) where N is a number of points

    :param initial_pressure: pressure in the reservoir before injection
    :type initial_pressure: numpy.array of shape (N, ) where N is a number of points

    :returns:
        out: Dictionary with keys 'pressure_front', 'plume_extent', and 'delineated_area'.

    """
    # Calculate difference in the initial and current pressures
    delta_pressure = pressure - initial_pressure

    # Calculate allowable pressure change
    allowable_delta_pressure = - initial_pressure + p['aquifer_pressure'] + \
        p['res_brine_density']*GRAVITY_CONSTANT*(p['reservoir_depth']-p['aquifer_depth'])

    # Calculate pressure front
    pressure_front = np.zeros(pressure.shape)
    pressure_front[np.where(delta_pressure > allowable_delta_pressure)] = 1

    # Calculate CO2 plume extent
    plume_extent = np.zeros(saturation.shape)
    plume_extent[np.where(saturation > p['saturation_threshold'])] = 1

    delineated_area = np.logical_or(pressure_front, plume_extent)

    # Setup dictionary of outputs
    out = {'pressure_front': pressure_front,
           'plume_extent': plume_extent,
           'delineated_area': delineated_area}

    return out


def hydrostatic_case_delineation(p, pressure, saturation, initial_pressure):
    """
    Calculate pressure front, CO2 plume extent and delineated area

    :param p: input parameters
    :type p: dict()

    :param pressure: pressure in the reservoir
    :type pressure: numpy.array of shape (N, ) where N is a number of points

    :param saturation: CO2 saturation in the reservoir
    :type saturation: numpy.array of shape (N, ) where N is a number of points

    :param initial_pressure: pressure in the reservoir before injection
    :type initial_pressure: numpy.array of shape (N, ) where N is a number of points

    :returns:
        out: Dictionary with keys 'pressure_front', 'plume_extent', and 'delineated_area'.

    """
    # Calculate difference in the initial and current pressures
    delta_pressure = pressure - initial_pressure

    # Calculate allowable pressure change
    allowable_delta_pressure = 0.5*GRAVITY_CONSTANT*(
        p['res_brine_density'] - p['aq_brine_density'])*(
            p['reservoir_depth'] - p['aquifer_depth'])

    # Calculate pressure front
    pressure_front = np.zeros(pressure.shape)
    pressure_front[np.where(delta_pressure > allowable_delta_pressure)] = 1

    # Calculate CO2 plume extent
    plume_extent = np.zeros(saturation.shape)
    plume_extent[np.where(saturation > p['saturation_threshold'])] = 1

    delineated_area = np.logical_or(pressure_front, plume_extent)

    # Setup dictionary of outputs
    out = {'pressure_front': pressure_front,
           'plume_extent': plume_extent,
           'delineated_area': delineated_area}

    return out


DELINEATION_FUNCTION = {1: simple_delineation,
                        2: under_pressurized_case_delineation,
                        3: hydrostatic_case_delineation}

class AreaEstimate(iam_bc.ComponentModel):
    """
    NRAP-Open-IAM Area Estimate component class.

    Area estimate component estimates 3 metrics applicable to the area of review:
    pressure front, CO2 plume extent and delinated area combining the first two.
    Depending on the initial conditions in the reservoir (e.g., underpressured
    or hydrostatic) different approaches to thresholds to the change in pressure
    can be applied.
    """
    def __init__(self, name, parent, coordinates=None, cell_size=1,
                 criteria=1, data_dim=2):
        """
        Constructor method of AreaEstimate class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param coordinates: coordinates associated with data
        :type coordinates: a numpy.array of shape (N, 2) for 2d data and
            np.array of shape (N, 3) for 3d data.

        :param cell_size: area or volume associated with each point/cell of data.
            It is used to estimate the size of plume in area or volume terms.
        :type cell_size: float or a numpy.array of length N

        :param criteria: flag variable indicating how the areas will be estimated
            Possible values:
                1 - simplest case; data will be compared to a provided threshold
                2 - pressure threshold will be calculated based on provided
                input parameters; applicable to under-pressurized case
                3 - pressure threshold will be calculated based on provided
                input parameters; applicable to hydrostatic initial case
        :type criteria: int

        :param data_dim: number of dimensions of data. Possible values:
            2 for 2d data and 3 for 3d data. Can also be obtained from
            coordinates if provided
        :type data_dim: int

        :returns: AreaEstimate class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        self.grid_obs_keys = AE_OBSERVATIONS + [f'max_{obs_nm}' for obs_nm in AE_OBSERVATIONS]

        # Add type attribute
        self.class_type = 'AreaEstimate'

        # Save attributes
        self.coordinates = coordinates
        self.cell_size = cell_size
        self.process_criteria(criteria)

        # Set default parameters of the component model
        # The next 2 parameters are used if criteria == 1
        self.add_default_par('delta_pressure_threshold', value=1.0e+6)
        self.add_default_par('saturation_threshold', value=1.0e-10)
        # The next parameters are used if criteria == 2
        # Parameter saturation_threshold is also used
        self.add_default_par('res_brine_density', value=1000.0)
        self.add_default_par('aquifer_depth', value=800.0)
        self.add_default_par('aquifer_pressure', value=5.0e+5)
        self.add_default_par('reservoir_depth', value=1000.0)
        # The next parameters are used if criteria == 3
        # Parameters res_brine_density, aquifer_depth, and reservoir_depth
        # are also used
        self.add_default_par('aq_brine_density', value=980.0)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['delta_pressure_threshold'] = [0, np.inf]
        self.pars_bounds['saturation_threshold'] = [0.0, 1.0]
        self.pars_bounds['res_brine_density'] = [800, 1100]
        self.pars_bounds['aq_brine_density'] = [800, 1100]
        self.pars_bounds['aquifer_depth'] = [1, 1500]
        self.pars_bounds['reservoir_depth'] = [500, 10000]
        self.pars_bounds['aquifer_pressure'] = [1.0e+5, 6.0e+7]

        # Add accumulators to keep track of intermediate plume locations
        for obs_nm in AE_OBSERVATIONS:
            self.add_accumulator(obs_nm, sim=0.0)
            self.add_accumulator('max_'+obs_nm, sim=0.0)

        debug_msg = 'AreaEstimate component created with name {}'.format(name)
        logging.debug(debug_msg)

    def process_criteria(self, value):
        if value in [1, 2, 3]:
            self.criteria = value
        else:
            err_msg = ''.join([
                "Argument 'criteria' of Area Estimate constructor method ",
                "has a wrong value of {}"]).format(value)
            raise ValueError(err_msg)

    def check_data_size(self, data):
        if self.coordinates is not None:
            for key in data:
                if len(data[key]) != len(self.coordinates):
                    warn_msg = "".join(["Number of data points provided through ",
                                        "argument 'coordinates' does not equal ",
                                        "to the size of '{}' data."]).format(key)
                    logging.warning(warn_msg)

    def update_all_parameters(self, scalar_p, **kwarg_p):
        """
        Transform all model parameters into array versions of themselves.
        """
        # Update major arrays
        combined_p = {}
        for key in PARAMETERS_GROUPS[self.criteria]:
            if key in kwarg_p:
                if len(kwarg_p[key]) != self.num_cells:
                    err_msg = "".join(["Length of input argument '{}' ({}) does not ",
                                       "equal to the number of data points {}"]).format(
                                           key, len(kwarg_p[key]), self.num_cells)
                    logging.error(err_msg)
                    raise IndexError(err_msg)
                combined_p[key] = np.array(kwarg_p[key])
            elif key in scalar_p:
                combined_p[key] = scalar_p[key]*np.ones(self.num_cells)

        return combined_p

    def simulation_model(self, p, time_point=365.25,
                         pressure=None, CO2saturation=None,
                         initial_pressure=None, **kwargs):
        """
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this module's docstring.

        :param p: input parameters of the component
        :type p: dict

        :param pressure: pressure in the injection zone/reservoir, in Pa;
            by default, its value is None
        :type pressure: numpy.array of floats of shape (N,) where N is number of
            data points

        :param CO2saturation: saturation of |CO2| phase in the injection zone;
            by default, its value is None
        :type CO2saturation: numpy.array of floats of shape (N,) where N is number of
            data points

        :param initial_pressure: pressure in the injection zone/reservoir
            before the injection (at time 0), in Pa; by default, its value is None
        :type initial_pressure: numpy.array of floats of shape (N,) where N is number of
            data points

        :param time_point: time point in days at which the estimates are to be
            made
        :type time_point: float

        :param kwargs: additional keyword arguments of the method. Possible keys
            are 'aquifer_pressure', 'aquifer_depth', 'aq_brine_density',
            'reservoir_depth', 'res_brine_density'

            'aquifer_pressure': pressure in the USDW zone; used when criteria
            is 2; numpy.array of floats of shape (N,) where N is number of
            data points or scalar if applicable

            'aquifer_depth': (representative) depth of the aquifer/USDW; used
            when criteria is 2 or 3; numpy.array of floats of shape (N,) where
            N is number of data points, or scalar if applicable

            'aq_brine_density': brine density in the aquifer/USDW

            'reservoir_depth': (representative) depth of the reservoir/injection
            zone; used when criteria is 2 or 3; numpy.array of floats of
            shape (N,) where N is number of data points, or scalar if applicable

            'res_brine_density': brine density in the injection zone

        :returns: out - dictionary of observations of area estimate component;
            keys: ['pressure_front', 'plume_extent', 'delineated_area']
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Check data size for pressure/saturation vs coordinates
        self.check_data_size({'pressure': pressure,
                              'CO2saturation': CO2saturation,
                              'initial_pressure': initial_pressure})

        self.num_cells = len(pressure)
        combined_p = self.update_all_parameters(actual_p, **kwargs)

        # Check how area has to be delineated
        if self.criteria in [1, 2, 3]:
            params = {k: combined_p[k] for k in PARAMETERS_GROUPS[self.criteria]}
            out = DELINEATION_FUNCTION[self.criteria](
                params, pressure, CO2saturation, initial_pressure)

        # Update accumulators related to the current time point
        for obs_nm in AE_OBSERVATIONS:
            self.accumulators[obs_nm].sim = out[obs_nm]

        # Update accumulators
        time_index = np.where(self._parent.time_array==time_point)[0][0]
        if time_index == 0:
            for obs_nm in AE_OBSERVATIONS:
                self.accumulators['max_'+obs_nm].sim = out[obs_nm]
                out['max_'+obs_nm] = out[obs_nm]
        else:
            for obs_nm in AE_OBSERVATIONS:
                self.accumulators['max_'+obs_nm].sim = np.logical_or(
                    self.accumulators['max_'+obs_nm].sim, out[obs_nm])
                out['max_'+obs_nm] = self.accumulators['max_'+obs_nm].sim

        out['area'] = np.sum(out['delineated_area']*self.cell_size)
        out['max_area'] = np.sum(out['max_delineated_area']*self.cell_size)

        # Return dictionary of outputs
        return out

    def reset(self):
        return
