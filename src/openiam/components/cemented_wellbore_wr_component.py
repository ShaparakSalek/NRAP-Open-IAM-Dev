# -*- coding: utf-8 -*-
import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt

from openiam.matk.ordereddict import OrderedDict

try:
    import openiam.components.iam_base_classes as iam_bc
except ImportError as err:
    print('Unable to load NRAP-Open-IAM base classes module: {}'.format(err))

from openiam.cf_interface.commons import process_parameters, process_dynamic_inputs
from openiam.cf_interface.strata import (get_comp_types_strata_pars,
                                         get_comp_types_strata_obs,
                                         get_strat_param_dict_for_link)

try:
    import openiam.components.models.wellbore.cemented as cwmodel
    import openiam.components.models.wellbore.cemented.cemented_wellbore_wr_ROM as wrrom
    from openiam.components.models.wellbore.cemented.utilities import (
        first_derivative, second_derivative)
except ImportError:
    print('\nERROR: Unable to load ROM for Cemented Wellbore WR component\n')
    sys.exit()


class CementedWellboreWR(iam_bc.ComponentModel):
    """
    The Cemented Wellbore WR (wider ranges) component model is an updated version
    of cemented wellbore component available in NRAP-Open-IAM. The model is
    built off detailed full-physics Finite Element Heat and Mass (FEHM)
    simulations, :cite:`Zyvoloski2007`. The FEHM simulations are
    three-dimensional (3-D), multiphase solutions of heat and mass transfer of
    water and supercritical, liquid, and gas |CO2|. After the simulations are
    completed, the surrogate model is built based on the key input parameters
    and corresponding output parameters. The approximate (surrogate) model is
    represented by polynomials in terms of input parameters that then can be
    sampled to estimate leakage rate for wells.  Early development work can be found
    in :cite:`RN1606`.

    When using the control file interface with more than 3 shale layers,
    the ``ThiefZone`` keyword can be used to specify the thief zone aquifer and the
    ``LeakTo`` keyword can be specified to name the upper aquifer.  These values will
    default to *aquifer1* and *aquifer2* respectively.
    In the FEHM simulations used to create the surrogate model some of the stratigraphy
    layers were setup with a fixed thickness. In particular, shale above aquifer
    had thickness 9.6 m; shale layer between the aquifer and the thief zone varied
    between 228.4 m and 2902.9 m.

    Component model input definitions:

    * **logWellPerm** [|log10| |m^2|] (-13.95 to -10.1) - logarithm of wellbore
      permeability (default: -13)

    * **logThiefPerm** [|log10| |m^2|] (-13.986 to -12.023) - logarithm of thief zone
      permeability (default: -12.5)

    * **thiefZoneThickness** [|m|] (0 to 99) - thickness of thief zone (default: 50);
      in 2 aquifers system linked to Stratigraphy parameter aquifer1Thickness

    * **aquiferThickness** [|m|] (10 to 230) - thickness of aquifer (default: 50);
      in 2 aquifers system linked to Stratigraphy parameter aquifer2Thickness

    * **reservoirThickness** [|m|] (16 to 193) - thickness of reservoir (default: 50);
      *linked to Stratigraphy*

    * **wellRadius** [|m|] (0.025 to 0.25) - radius of the wellbore (default: 0.05)

    * **initPressure** [|Pa|] (1.0e+5 to 5.0e+7) - initial pressure at the base of
      the wellbore (default: 2.0e+7 |Pa|, or 20 |MPa|); *from linked component*

    * **wellDepth** [|m|] (1056.0 to 3194.0) - depth in meters from ground surface to
      top of reservoir (default: 1500); *linked to Stratigraphy*

    * **depthRatio** [-] (0.054 to 0.69798) - fraction of well depth to the center of
      the thief zone from the top of the reservoir (default: 0.5);
      *linked to Stratigraphy*.

    Temporal inputs of the Cemented Wellbore component are not provided directly
    to the component model method but rather are calculated from the current
    and several past values of pressure and |CO2| saturation.
    The calculated temporal inputs are then checked against the boundary assumptions
    of the underlying reduced order model. The Cemented Wellbore component model
    temporal inputs are:

    * **deltaP** [|Pa|] (-1246850.9 to 8326262.6) - difference between the current
      and initial pressure at the wellbore

    * **pressurePrime** [|Pa/s|] (-5665.23 to 8898.32) - first pressure derivative

    * **pressureDPrime** [|Pa/s^2|] (-232.64 to 232.58) - second pressure derivative

    * **saturation** [-] (4.63e-10 to 1.0) - |CO2| saturation at the wellbore

    * **saturationPrime** [|1/s|] (-2.638e-7 to 1.294e-3) - first |CO2|
      saturation derivative

    * **saturationDPrime** [|1/s^2|] (-8.806e-6 to 7.269e-6) - second |CO2|
      saturation derivative.

    The possible outputs from the Cemented Wellbore component are leakage rates
    of |CO2| and brine to aquifer, thief zone and atmosphere. The names of the
    observations are of the form:

    * **CO2_aquifer1**, **CO2_aquifer2**, **CO2_atm** [|kg/s|] - |CO2| leakage rates

    * **brine_aquifer1**, **brine_aquifer2**, **brine_atm** [|kg/s|] -
      brine leakage rates

    * **mass_CO2_aquifer1**, **mass_CO2_aquifer2** [|kg|] -
      mass of |CO2| leaked into aquifers.

    """
    def __init__(self, name, parent, header_file_dir=None, ignore_tmp_input_check=True):
        """
        Constructor method of CementedWellbore class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param header_file_dir: path to the directory that contains
            the lookup tables
        :type workdir: str

        :param ignore_tmp_inp_check: flag indicating whether a user wants to ignore
            the messages about the temporal inputs not satisfying the boundaries of the ROM,
            allows truncation and wants to continue running the component. Otherwise,
            error message will be shown and simulation will stop.

        :returns: CementedWellboreWR class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'CementedWellboreWR'

        # Define output dictionary labels
        self.output_labels = ['CO2_aquifer1', 'CO2_aquifer2', 'CO2_atm',
                              'brine_aquifer1', 'brine_aquifer2', 'brine_atm',
                              'mass_CO2_aquifer1', 'mass_CO2_aquifer2']

        # Setup default observations of the component
        self.default_obs = {obs_nm: 0.0 for obs_nm in self.output_labels}

        # Set default parameters of the component model
        self.add_default_par('wellDepth', value=1500.0)
        self.add_default_par('depthRatio', value=0.5)
        self.add_default_par('aquiferThickness', value=50.0)
        self.add_default_par('thiefZoneThickness', value=50.0)
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('logWellPerm', value=-13.0)
        self.add_default_par('logThiefPerm', value=-12.5)
        self.add_default_par('initPressure', value=2.0e+7)
        self.add_default_par('wellRadius', value=0.05)
        self.add_default_par('g', value=9.8)

        # Check whether header_file_dir is provided
        if header_file_dir is None:
            header_file_dir = os.sep.join([cwmodel.__path__[0], 'header_files_wr'])

        # Instantiate solution object
        self.sol = wrrom.Solution(header_file_dir)

        # Update boundaries of the log well permeability parameter
        for nm in self.sol.roms:
            self.sol.roms[nm].change_boundaries(2, -13.95, -10.1)

        rom = self.sol.roms['CO2ThiefZone']  # pick any ROM
        # Define dictionary of model parameters boundaries
        # Boundaries of some parameters are defined by the ROM
        self.pars_bounds = dict()
        self.pars_bounds['wellDepth'] = [rom._mins[0], rom._maxs[0]]     #  [1056.0, 3194.0]
        self.pars_bounds['depthRatio'] = [rom._mins[1], rom._maxs[1]]    #  [0.05403158769742311, 0.6979895104895105]
        self.pars_bounds['logWellPerm'] = [rom._mins[2], rom._maxs[2]]   #  [-13.95, -10.1]
        self.pars_bounds['logThiefPerm'] = [rom._mins[3], rom._maxs[3]]  #  [-13.986, -12.023]
        self.pars_bounds['thiefZoneThickness'] = [rom._mins[4], rom._maxs[4]] # [0.0, 99.0]
        self.pars_bounds['aquiferThickness'] = [rom._mins[5], rom._maxs[5]]   # [10.0, 230.0]
        self.pars_bounds['reservoirThickness'] = [rom._mins[6], rom._maxs[6]] # [16.0, 193.0]
        self.pars_bounds['initPressure'] = [1.0e+5, 5.0e+7]
        self.pars_bounds['wellRadius'] = [0.025, 0.25]
        self.pars_bounds['g'] = [9.76, 9.83]

        debug_msg = "Component {} parameters boundaries: {}".format(
            self.name, self.pars_bounds)
        logging.debug(debug_msg)

        # Define dictionary of temporal data limits
        # Boundaries of temporal inputs are defined by the ROM
        self.temp_data_bounds = dict()
        MPa = 1.0e+6
        self.temp_data_bounds['deltaP'] = [
            'change in pressure',
            rom._mins[7]*MPa, rom._maxs[7]*MPa] #  [-1246850.999999999, 8326262.600000003] Pa

        self.temp_data_bounds['pressurePrime'] = [
            'first pressure derivative',
            rom._mins[8]*MPa, rom._maxs[8]*MPa]  #  [-5665.231670491607, 8898.320328542035] Pa/s

        self.temp_data_bounds['pressureDPrime'] = [
            'second pressure derivative',
            rom._mins[9]*MPa, rom._maxs[9]*MPa]  #  [-232.6444260521089, 232.58222171184337] Pa/s^2

        self.temp_data_bounds['saturation'] = [
            'saturation',
            rom._mins[10], rom._maxs[10]]  #  [4.6354276e-10, 1.0]

        self.temp_data_bounds['saturationPrime'] = [
            'first saturation derivative',
            rom._mins[11], rom._maxs[11]]  #  [-2.628459958932679e-07, 0.001294187515400389] 1/s

        self.temp_data_bounds['saturationDPrime'] = [
            'second saturation derivative',
            rom._mins[12], rom._maxs[12]]  #  [-8.806248708726336e-06, 7.268626818201906e-06] 1/s^2

        debug_msg = "Component {} temporal inputs boundaries: {}".format(
            self.name, self.temp_data_bounds)
        logging.debug(debug_msg)

        # Set the flag reflecting user attitude toward temporal input check
        self.ignore_tmp_input_check = ignore_tmp_input_check

        # Add accumulators for CO2 mass, time, saturation and pressure
        for i in range(2):
            self.add_accumulator('mass_CO2_aquifer'+str(i+1), sim=0.0)
            # time_point_nmin is a previous time point (in days);
            # by default, it's -1.0 to indicate that pressure/saturation data
            # for this point are not available
            self.add_accumulator('time_point_nmin'+str(i+1), sim=-1.0)
            self.add_accumulator('pressure_nmin'+str(i+1), sim=0.0)
            self.add_accumulator('CO2saturation_nmin'+str(i+1), sim=0.0)

        # Default indices for aquifer and thief zone
        self.leak_layer = 2
        self.thief_zone = 1
        self.default_out = {}

        self.init_pressure = False

        debug_msg = 'CementedWellboreWR component created with name {}'.format(name)
        logging.debug(debug_msg)

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        for key, val in p.items():
            if key in self.pars_bounds:
                if (val < self.pars_bounds[key][0]) or (
                        val > self.pars_bounds[key][1]):
                    if key in ['wellDepth', 'depthRatio']:
                        error_msg = 'Current {} parameter value: {}'.format(key, val)
                        logging.error(error_msg)
                        error_msg = ''.join([
                            'Parameter {} of CementedWellboreWR is out of boundaries. ',
                            'The reservoir and/or well configuration do not satisfy ',
                            'assumptions of component.']).format(key)
                        logging.error(error_msg)
                        raise ValueError(error_msg)

                    warn_msg = ''.join([
                        'Parameter {} of CementedWellboreWR component {} ',
                        'is out of boundaries.']).format(key, self.name)
                    error_msg = 'Current {} parameter value: {} and boundaries are [{}, {}].'.format(
                        key, val, self.pars_bounds[key][0], self.pars_bounds[key][1])
                    logging.error(error_msg)
                    logging.warning(warn_msg)
                    if key in ['logWellPerm', 'logThiefPerm']:
                        warn_msg = ''.join([
                            'Parameter {} of CementedWellboreWR component {} ',
                            'will be modified to satisfy the component ',
                            'requirements.']).format(key, self.name)
                        logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of CementedWellboreWR component {}.']).format(key, self.name)
                logging.warning(warn_msg)

    def check_temporal_inputs(self, time, temp_inputs):
        """
        Check whether temporal data fall within specified boundaries.

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict
        """
        debug_msg = 'Temporal inputs of component {} are {}'.format(
            self.name, temp_inputs)
        logging.debug(debug_msg)

        for key, val in temp_inputs.items():
            if (val < self.temp_data_bounds[key][1]) or (
                    val > self.temp_data_bounds[key][2]):
                if self.ignore_tmp_input_check:
                    warn_msg1 = ''.join([
                        'Temporal input {} ({}) of CementedWellboreWR component {} ',
                        'is outside the model range [{}, {}] at time t = {} days']).format(
                            self.temp_data_bounds[key][0].lower(), val, self.name,
                            self.temp_data_bounds[key][1],
                            self.temp_data_bounds[key][2], time)
                    logging.warning(warn_msg1)
                    warn_msg2 = ''.join([
                        'Temporal input {} will be modified to satisfy ',
                        'the ROM boundaries.']).format(
                            self.temp_data_bounds[key][0].lower())
                    logging.warning(warn_msg2)
                else:
                    warn_msg1 = ''.join([
                        'Temporal input {} ({}) of component {} is outside ',
                        'the model range at time t = {} days']).format(
                            self.temp_data_bounds[key][0].lower(), val, self.name, time)
                    logging.warning(warn_msg1)
                    warn_msg2 = ''.join([
                        'User can change temporal input check preferences ',
                        'of CementedWellboreWR component by setting the component ',
                        'attribute ignore_tmp_input_check to True.'])
                    logging.error(warn_msg2)
                    raise ValueError('Temporal inputs are outside the ROM limits.')

    def simulation_model(self, p, time_point=365.25, time_step=365.25,
                         pressure=0.0, CO2saturation=0.0):
        """
        Return CO2 and brine leakage rates corresponding to the provided input.
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this modules docstring.

        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: current time point in days (at which the output is
            to be provided); by default its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between current and previous time points
            in days; by default its value is 365.25 (1 year in days)
        :type time_point: float

        :param pressure: pressure at the bottom of leaking well
            at the current time point, in Pa; by default its value is 0.0
        :type pressure: float

        :param saturation: CO2 saturation at the bottom of leaking well;
            by default its value is 0.0
        :type saturation: float

        :returns: dictionary of observations (leakage rates, etc.) of cemented
            wellbore model; keys:
            ['CO2_aquifer1','CO2_aquifer2','CO2_atm',
            'brine_aquifer1','brine_aquifer2','brine_atm',
            ''mass_CO2_aquifer1','mass_CO2_aquifer2']
        """
        # Return the initial state of the model if requested
        if time_point == 0.0:
            # For time point 0 all outputs of the model are zeros
            self.init_pressure = pressure
            out = self.default_out.copy()
            out.update(dict(list(zip(self.output_labels, len(self.output_labels)*[0.0]))))
            # Setup accumulators for the next time point
            self.accumulators['time_point_nmin1'].sim = 0.0
            self.accumulators['time_point_nmin2'].sim = -1.0
            self.accumulators['pressure_nmin1'].sim = pressure
            self.accumulators['pressure_nmin2'].sim = 0.0
            self.accumulators['CO2saturation_nmin1'].sim = CO2saturation
            self.accumulators['CO2saturation_nmin2'].sim = 0.0
            self.accumulators['mass_CO2_aquifer1'].sim = 0.0
            self.accumulators['mass_CO2_aquifer2'].sim = 0.0
            return out

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Set parameters
        wellDepth = actual_p['wellDepth']
        depthRatio = actual_p['depthRatio']
        logWellPermeability = actual_p['logWellPerm']
        logThiefZonePermeability = actual_p['logThiefPerm']
        aquiferThickness = actual_p['aquiferThickness']
        thiefZoneThickness = actual_p['thiefZoneThickness']
        reservoirThickness = actual_p['reservoirThickness']
        if self.init_pressure:
            initialPressure = self.init_pressure # actual_p['initPressure']
        else:
            initialPressure = actual_p['initPressure']

        # Get values from accumulator variables
        time_point_nmin1 = self.accumulators['time_point_nmin1'].sim
        time_point_nmin2 = self.accumulators['time_point_nmin2'].sim
        pressure_nmin1 = self.accumulators['pressure_nmin1'].sim
        pressure_nmin2 = self.accumulators['pressure_nmin2'].sim
        CO2saturation_nmin1 = self.accumulators['CO2saturation_nmin1'].sim
        CO2saturation_nmin2 = self.accumulators['CO2saturation_nmin2'].sim
        CO2MassInAquifer1 = self.accumulators['mass_CO2_aquifer1'].sim
        CO2MassInAquifer2 = self.accumulators['mass_CO2_aquifer2'].sim

        # Check which keyword arguments are provided
        if time_point_nmin1 > 0.0 and time_point_nmin2 >= 0.0:   # three data points are available
            deltaP = pressure-initialPressure
            pressurePrime = first_derivative(pressure, pressure_nmin1, time_step)
            pressureDPrime = second_derivative(
                pressure, pressure_nmin1, pressure_nmin2,
                time_step, time_point_nmin1-time_point_nmin2)
            saturationPrime = first_derivative(
                CO2saturation, CO2saturation_nmin1, time_step)
            saturationDPrime = second_derivative(
                CO2saturation, CO2saturation_nmin1, CO2saturation_nmin2,
                time_step, time_point_nmin1-time_point_nmin2)
            saturation = CO2saturation
        elif time_point_nmin1 >= 0.0:    # the current and previous points are available
            deltaP = pressure-initialPressure
            pressurePrime = first_derivative(pressure, pressure_nmin1, time_step)
            pressureDPrime = 0.0
            saturationPrime = first_derivative(
                CO2saturation, CO2saturation_nmin1, time_step)
            saturationDPrime = 0.0
            saturation = CO2saturation

        # Define array of input parameters
        megaPascal = 1.0e+6
        sec_per_day = 24*3600
        scalar1 = megaPascal*sec_per_day
        scalar2 = megaPascal*sec_per_day**2

        # Check whether pressure and saturation inputs satisfy the model requirements
        self.check_temporal_inputs(time_point, dict(list(zip(
            ['deltaP', 'pressurePrime', 'pressureDPrime',
             'saturation', 'saturationPrime', 'saturationDPrime'],
            [deltaP, pressurePrime/sec_per_day, pressureDPrime/sec_per_day**2,
             saturation, saturationPrime/sec_per_day, saturationDPrime/(sec_per_day**2)]))))

        inputArray = np.array([
            wellDepth, depthRatio, logWellPermeability, logThiefZonePermeability,
            thiefZoneThickness, aquiferThickness, reservoirThickness,
            deltaP/megaPascal, pressurePrime/scalar1, pressureDPrime/scalar2,
            saturation, saturationPrime/sec_per_day, saturationDPrime/(sec_per_day**2)])
        self.sol.find(inputArray)

        # Rescale solution by constant dependent on the wellbore radius
        # TODO figure out the radius of the well used in the simulation
        wellRadius = actual_p['wellRadius']
        originalRadius = 0.05
        solScalar = (wellRadius/originalRadius)**2
        CO2LeakageRates = solScalar*self.sol.CO2LeakageRates
        brineLeakageRates = solScalar*self.sol.brineLeakageRates

        # Set values to leakage rates
        out = self.default_out.copy()
        out.update(dict(list(zip(self.output_labels[0:6], np.concatenate(
            [CO2LeakageRates, brineLeakageRates])))))
        self.accumulators['time_point_nmin1'].sim = time_point
        self.accumulators['time_point_nmin2'].sim = time_point_nmin1
        self.accumulators['pressure_nmin1'].sim = pressure
        self.accumulators['pressure_nmin2'].sim = pressure_nmin1
        self.accumulators['CO2saturation_nmin1'].sim = CO2saturation
        self.accumulators['CO2saturation_nmin2'].sim = CO2saturation_nmin1

        out['mass_CO2_aquifer{tz}'.format(tz=self.thief_zone)] = max(
            0.0, CO2MassInAquifer1 + (
                CO2LeakageRates[0]-CO2LeakageRates[1])*time_step*sec_per_day)
        out['mass_CO2_aquifer{ll}'.format(ll=self.leak_layer)] = max(
            0.0, CO2MassInAquifer2 + (
                CO2LeakageRates[1]-CO2LeakageRates[2])*time_step*sec_per_day)
        self.accumulators['mass_CO2_aquifer1'].sim = out['mass_CO2_aquifer{tz}'.format(
            tz=self.thief_zone)]
        self.accumulators['mass_CO2_aquifer2'].sim = out['mass_CO2_aquifer{ll}'.format(
            ll=self.leak_layer)]

        return out

    def connect_with_system(self, component_data, name2obj_dict, *args, **kwargs):
        """
        Code to add cemented wellbore to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :returns: None
        """
        # Process parameters data if any
        process_parameters(self, component_data)

        # Process dynamic inputs if any
        process_dynamic_inputs(self, component_data)

        # These lists indicate the stratigraphy component types that offer thicknesses
        # and depths as parameters or as observations.
        types_strata_pars = get_comp_types_strata_pars()
        types_strata_obs = get_comp_types_strata_obs()

        # Get strata component
        strata = name2obj_dict['strata']

        # Determine number of shale layers
        num_shale_layers = strata.get_num_shale_layers()

        # Determine whether it's an unexpected scenario with more than 2 aquifers
        if 'LeakTo' in component_data:
            leak_to = component_data['LeakTo'].lower()
        else:
            leak_to = 'aquifer2'

        if 'aquifer' in leak_to:
            self.leak_layer = int(leak_to[7:])
        else:
            self.leak_layer = num_shale_layers - 1

        # Determine to which aquifer the thief zone is assigned to
        if 'ThiefZone' in component_data:
            thief_zone = component_data['ThiefZone'].lower()
        else:
            thief_zone = 'aquifer1'

        self.thief_zone = int(thief_zone[7:])

        # wellDepth can also come from table, so we need to add check for the
        # presence of the parameter among parameters of the component itself
        if name2obj_dict['strata_type'] in types_strata_pars:
            if 'wellDepth' not in self.pars and \
                    'wellDepth' not in self.deterministic_pars:
                res_depth = '{}.shale1Thickness '.format(strata.name)
                for il in range(1, num_shale_layers):
                    res_depth += ('+ {nm}.aquifer{il}Thickness + {nm}.shale{ilp1}Thickness '
                                  .format(nm=strata.name, il=il, ilp1=il+1))
                self.add_composite_par('wellDepth', res_depth)

        elif name2obj_dict['strata_type'] in types_strata_obs:
            if 'wellDepth' in self.pars or 'wellDepth' in self.deterministic_pars:
                warning_msg = strata.parameter_assignment_warning_msg().format(
                    'wellDepth', 'shale1Depth')

                logging.warning(warning_msg)

                if 'wellDepth' in self.pars:
                    del self.pars['wellDepth']

                elif 'wellDepth' in self.deterministic_pars:
                    del self.deterministic_pars['wellDepth']

            self.add_par_linked_to_obs('wellDepth', strata.linkobs['shale1Depth'])

        tz = self.thief_zone
        if 'depthRatio' in self.pars or 'depthRatio' in self.deterministic_pars:
            warn_msg = ''.join(['Parameter depthRatio is defined by user and, ',
                                'thus, is not linked to stratigraphy. ',
                                'In general, it is not a recommended setting. ',
                                'Please check your setup.'])
            logging.warning(warn_msg)
        else:
            depth_ratio = '({nm}.aquifer{tz}Thickness/2'.format(
                nm=strata.name, tz=tz)
            for il in range(tz+1, num_shale_layers):
                depth_ratio += ' + {nm}.shale{il}Thickness + {nm}.aquifer{il}Thickness'.format(
                    nm=strata.name, il=il)
            depth_ratio += ' + {nm}.shale{nsl}Thickness)'.format(
                nm=strata.name, nsl=num_shale_layers)
            depth_ratio += '/{selfm}.wellDepth'.format(selfm=self.name)

            if name2obj_dict['strata_type'] in types_strata_pars:
                self.add_composite_par('depthRatio', depth_ratio)

            elif name2obj_dict['strata_type'] in types_strata_obs:
                self.add_par_linked_to_composite_obs('depthRatio', depth_ratio)

        # Setup the rest of the stratigraphy related parameters
        setup_data = {'aquifer': ['aquifer', self.leak_layer],
                      'thiefZone': ['aquifer', tz],
                      'reservoir': ['reservoir', '']}

        for key, vals in setup_data.items():
            par_name = '{}Thickness'.format(key)
            if par_name not in self.pars and par_name not in self.deterministic_pars:
                strata_par_name = '{}{}Thickness'.format(vals[0], vals[1])

                if name2obj_dict['strata_type'] in types_strata_pars:
                    connect = get_strat_param_dict_for_link(strata_par_name, strata)

                    if connect is None:
                        strata_par_name = '{}Thickness'.format(vals[0])

                        connect = get_strat_param_dict_for_link(strata_par_name, strata)

                        if connect is None:
                            err_msg = ''.join([
                                'Unable to find "{}" or "{}" parameters. Please ',
                                'check setup of the cemented wellbore component ',
                                'and stratigraphy.']).format(par_name, strata_par_name)
                            logging.error(err_msg)
                            raise KeyError(err_msg)

                    self.add_par_linked_to_par(par_name, connect[strata_par_name])

                elif name2obj_dict['strata_type'] in types_strata_obs:
                    self.add_par_linked_to_obs(par_name, strata.linkobs[strata_par_name])

        # Make model connections
        if 'Connection' in component_data:
            connection = None
            try:
                connection = name2obj_dict[component_data['Connection']]
            except KeyError:
                pass
            system_inputs = ['pressure', 'CO2saturation']
            for sinput in system_inputs:
                connection.add_obs_to_be_linked(sinput)
                self.add_kwarg_linked_to_obs(sinput, connection.linkobs[sinput])

        for il in range(1, num_shale_layers):
            aq = 'aquifer{}'.format(il)
            self.default_out['CO2_' + aq] = 0.0
            self.default_out['brine_' + aq] = 0.0
        self.output_labels = ['CO2_aquifer{tz}', 'CO2_aquifer{ll}', 'CO2_atm',
                              'brine_aquifer{tz}', 'brine_aquifer{ll}', 'brine_atm',
                              'mass_CO2_aquifer{tz}', 'mass_CO2_aquifer{ll}']
        self.output_labels = [l.format(ll=self.leak_layer, tz=self.thief_zone)
                              for l in self.output_labels]

    # Attributes for system connections
    system_inputs = ['pressure',
                     'CO2saturation']
    composite_inputs = OrderedDict()
    composite_inputs['wellDepth'] = ('({strata}.shale1Thickness+{strata}.shale2Thickness+' +
                                     '{strata}.shale3Thickness+{strata}.aquifer1Thickness+' +
                                     '{strata}.aquifer2Thickness)')
    composite_inputs['depthRatio'] = (
        '({strata}.shale2Thickness+{strata}.shale3Thickness' +
        '+ {strata}.aquifer2Thickness + {strata}.aquifer1Thickness/2)/{selfm}.wellDepth')


def test_cemented_wellbore_wr_component():
    try:
        from openiam.components.analytical_reservoir_component import AnalyticalReservoir
    except ImportError as err:
        print('Unable to load NRAP-Open-IAM class module: {}'.format(err))

    logging.basicConfig(level=logging.WARNING)

    # Create system model
    num_years = 50.
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}     # time is given in days
    sm = iam_bc.SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('numberOfShaleLayers', value=3, vary=False)
    ares.add_par('shale1Thickness', min=500.0, max=550., value=525.0)
    ares.add_par('shale2Thickness', min=450.0, max=500., value=475.0)
    ares.add_par('shale3Thickness', value=35.0, vary=False)
    ares.add_par('aquifer1Thickness', value=45.0, vary=False)
    ares.add_par('aquifer2Thickness', value=23.0, vary=False)
    ares.add_par('reservoirThickness', value=55.0, vary=False)
    ares.add_par('injRate', value=0.1, vary=False)
    ares.add_par('reservoirRadius', value=1000, vary=False)

    # Add observations of reservoir component model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with name ares.obsnm_0, ares.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')

    # Add cemented wellbore component
    cw = sm.add_component_model_object(CementedWellboreWR(name='cw', parent=sm))

    # Add parameters of cemented wellbore component
    cw.add_par('logWellPerm', min=-13.9, max=-12.0, value=-13.0)
    cw.add_par_linked_to_par('thiefZoneThickness',
                             ares.deterministic_pars['aquifer1Thickness'])
    cw.add_par_linked_to_par('aquiferThickness',
                             ares.deterministic_pars['aquifer2Thickness'])
    cw.add_par_linked_to_par('reservoirThickness',
                             ares.deterministic_pars['reservoirThickness'])

    # Add keyword arguments of the cemented wellbore component model
    cw.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    cw.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])

    # Add composite parameters of cemented wellbore component
    #print(ares.pars['shale2Thickness'].name, ares.deterministic_pars['shale3Thickness'].name)
#    cw.add_composite_par('wellDepth', expr=ares.pars['shale1Thickness'].name+
#        '+'+ares.pars['shale2Thickness'].name+
#        '+'+ares.deterministic_pars['shale3Thickness'].name+
#        '+'+ares.deterministic_pars['aquifer1Thickness'].name+
#        '+'+ares.deterministic_pars['aquifer2Thickness'].name)
#    cw.add_composite_par('depthRatio',
#        expr='(' + ares.pars['shale2Thickness'].name +
#        '+' + ares.deterministic_pars['shale3Thickness'].name +
#        '+' + ares.deterministic_pars['aquifer2Thickness'].name +
#        '+' + ares.deterministic_pars['aquifer1Thickness'].name + '/2)/' +
#        cw.composite_pars['wellDepth'].name)
#    cw.add_composite_par('initPressure',
#        expr=ares.default_pars['datumPressure'].name +
#        '+' + cw.composite_pars['wellDepth'].name + '*98/10' +
#        '*' + ares.default_pars['brineDensity'].name)

    # Shorter way to write the same thing as the commented code above
    cw.add_composite_par('wellDepth', expr='ares.shale1Thickness' +
                         '+ares.shale2Thickness + ares.shale3Thickness' +
                         '+ares.aquifer1Thickness+ ares.aquifer2Thickness')
    cw.add_composite_par('depthRatio',
                         expr='(ares.shale2Thickness+ares.shale3Thickness' +
                         '+ ares.aquifer2Thickness + ares.aquifer1Thickness/2)/cw.wellDepth')
    cw.add_composite_par('initPressure',
                         expr='ares.datumPressure + cw.wellDepth*cw.g*ares.brineDensity')

    # Add observations of the cemented wellbore component
    cw.add_obs('CO2_aquifer1')
    cw.add_obs('CO2_aquifer2')
    cw.add_obs('CO2_atm')
    cw.add_obs('brine_aquifer1')
    cw.add_obs('brine_aquifer2')
    cw.add_obs('mass_CO2_aquifer2')

    # Show model parameters and temporal inputs boundaries
    cw.show_input_limits()

    # Run system model forward
    sm.forward()

    # Since the observations at the particular time points are different variables,
    # method collect_observations_as_time_series creates lists of
    # values of observations belonging to a given component (e.g. cw) and having the same
    # base name (e.g. 'CO2_aquifer1', etc) but differing in indices.
    # More details are given in the docstring and documentation to the method
    # collect_observations_as_time_series of SystemModel class.
    CO2leakrates_aq1 = sm.collect_observations_as_time_series(cw, 'CO2_aquifer1')
    CO2leakrates_aq2 = sm.collect_observations_as_time_series(cw, 'CO2_aquifer2')
    pressure = sm.collect_observations_as_time_series(ares, 'pressure')
    CO2saturation = sm.collect_observations_as_time_series(ares, 'CO2saturation')
    mass_CO2_aquifer2 = sm.collect_observations_as_time_series(cw, 'mass_CO2_aquifer2')

    # Print observations at all time points
    print('---------------------------------------------------')
    print('Pressure', pressure, sep='\n')
    print('---------------------------------------------------')
    print('CO2 saturation', CO2saturation, sep='\n')
    print('---------------------------------------------------')
    print('CO2 leakage rates to thief zone', CO2leakrates_aq1, sep='\n')
    print('---------------------------------------------------')
    print('Mass of CO2 accumulated in aquifer', mass_CO2_aquifer2, sep='\n')

    # Plot some results
    plt.figure(1)
    plt.plot(sm.time_array/365.25, CO2leakrates_aq1, color='#000055',
             linewidth=1, label="thief zone")
    plt.plot(sm.time_array/365.25, CO2leakrates_aq2, color='#FF0066',
             linewidth=1, label="aquifer")
    plt.legend()
    plt.xlabel('Time, t (years)')
    plt.ylabel('Leakage rates, q (kg/s)')
    plt.title(r'Leakage of CO$_2$: Thief Zone and Aquifer')
    plt.show()
