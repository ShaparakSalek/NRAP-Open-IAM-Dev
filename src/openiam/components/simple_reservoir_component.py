# -*- coding: utf-8 -*-
from math import sqrt
import sys
import os
import logging
from warnings import warn
import numpy as np
import matplotlib.pyplot as plt

try:
    import openiam.components.iam_base_classes as iam_bc
except ImportError as err:
    print('Unable to load NRAP-Open-IAM base classes module: {}'.format(err))

try:
    import openiam.components.models.reservoir.simple.simple_reservoir_ROM as sresrom
except ImportError:
    print('\nERROR: Unable to load ROM for Simple Reservoir component\n')
    sys.exit()


class SimpleReservoir(iam_bc.ComponentModel):
    """
    Note: Simple Reservoir component will be deprecated in the future versions
    of NRAP-Open-IAM. We recommend using Analytical Reservoir component instead.
    In this case updating the parameters might be needed.

    The Simple Reservoir component model is a semi-analytical model for the
    reservoir. It is focused on flow across relatively large distances
    and does not take into account discrete features of the flow paths such as
    fractures, cracks, etc. The model is based on work of Nordbotten et al.,
    :cite:`C2011`. Further reading can be found in :cite:`N2005a`, :cite:`C2009a`,
    :cite:`N2009`, :cite:`N2011a`.

    In the NRAP-Open-IAM control file, the type name for the Simple Reservoir component is
    ``SimpleReservoir``. The description of the component's parameters is
    provided below:

    * **logResPerm** [|log10| |m^2|] (-14 to -9) - logarithm of reservoir
      permeability (default: -12)

    * **reservoirPorosity** [-] (0.01 to 1) - porosity of reservoir (default: 0.3)

    * **brineDensity** [|kg/m^3|] (900 to 1500) - density of brine phase
      (default: 1000)

    * **CO2Density** [|kg/m^3|] (100 to 1500) - density of |CO2| phase
      (default: 479)

    * **brineViscosity** [|Pa*s|] (1.0e-4 to 5.0e-3) - viscosity of brine phase
      (default: 2.535e-3)

    * **CO2Viscosity** [|Pa*s|] (1.0e-6 to 1.0e-4)  - viscosity of |CO2| phase
      (default: 3.95e-5)

    * **brineResSaturation** [-] (0 to 0.7) - residual saturation of brine phase
      (default: 0.1)

    * **compressibility** [|Pa^-1|] (5.0e-11 to 1.0e-9) - compressibility of brine
      and |CO2| phases (assumed to be the same for both phases) (default: 5.1e-11)

    * **injRate** [|m^3/s|] (1.0e-3 to 10) - |CO2| injection rate (default: 0.1)

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3); *linked to Stratigraphy*. The shale units must be
      separated by an aquifer.

    * **shaleThickness** [|m|] (1 to 1600) - thickness of shale layers (default
      250); *linked to Stratigraphy*. Thickness of shale layer 1, for example,
      can be defined by **shale1Thickness**; otherwise, shale layers for which
      the thickness is not defined will be assigned a default thickness.

    * **aquiferThickness** [|m|] (1 to 1600) - thickness of aquifers (default: 100);
      *linked to Stratigraphy*. Thickness of aquifer 1, for example, can be defined
      by **aquifer1Thickness**; otherwise, aquifers for which the thickness
      is not defined will be assigned a default thickness.

    * **reservoirThickness** [|m|] (1 to 1600) - thickness of reservoir (default: 50);
      *linked to Stratigraphy*

    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the
      system (default: 101,325); *linked to Stratigraphy*.

    Possible observations from the Simple Reservoir component are:

    * **pressure** [|Pa|] - pressure at top of the reservoir
      at the user defined location(s)

    * **CO2saturation** [-] - |CO2| saturation at the top of the reservoir at the
      user defined location(s)

    * **mass_CO2_reservoir** [|kg|] - mass of the |CO2| in the reservoir.

    """
    def __init__(self, name, parent, injX=0., injY=0., locX=100., locY=100.):
        """
        Constructor method of SimpleReservoir class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param injX: x-coordinate of the injector location
        :type injX: float

        :param injY: x-coordinate of the injector location
        :type injY: float

        :param locX: x-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locX: float

        :param locY: y-coordinate of the location for pressure and CO2 saturation
            to be calculated

        :type locY: float

        :returns: SimpleReservoir class object
        """
        # Throw a deprecation warning on initialization of the class.
        # Using of stacklevel=2 will show a warning as being raised from
        # whatever is calling the __init__ method rather from the __init__ itself.
        warn(f'\n{self.__class__.__name__} class will be deprecated. Use AnalyticalReservoir instead.',
             DeprecationWarning, stacklevel=2)
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'SimpleReservoir'

        # Set default parameters of the component model
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('shaleThickness', value=250.0)
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('logResPerm', value=-12.0)
        self.add_default_par('reservoirPorosity', value=0.3)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('brineDensity', value=1000.0)
        self.add_default_par('CO2Density', value=479.0)
        self.add_default_par('brineViscosity', value=2.535e-3)
        self.add_default_par('CO2Viscosity', value=3.95e-5)
        self.add_default_par('brineResSaturation', value=0.1)
        self.add_default_par('compressibility', value=5.1e-11)
        self.add_default_par('injRate', value=0.1)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['shaleThickness'] = [1.0, 1600.0]
        self.pars_bounds['aquiferThickness'] = [1.0, 1600.0]
        self.pars_bounds['reservoirThickness'] = [1.0, 1600.0]
        self.pars_bounds['logResPerm'] = [-14.0, -9.0]
        self.pars_bounds['reservoirPorosity'] = [0.01, 1.0]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]
        self.pars_bounds['brineDensity'] = [900.0, 1500.0]
        self.pars_bounds['CO2Density'] = [100.0, 1500.0]
        self.pars_bounds['brineViscosity'] = [1.0e-4, 5.0e-3]
        self.pars_bounds['CO2Viscosity'] = [1.0e-6, 1.0e-4]
        self.pars_bounds['brineResSaturation'] = [0.0, 0.7]
        self.pars_bounds['compressibility'] = [5.0e-11, 1.0e-9]
        self.pars_bounds['injRate'] = [1.0e-3, 10.0]

        # Specify default locations of injection and leaking wells
        self.injX = injX
        self.injY = injY
        self.locX = locX
        self.locY = locY

        # Setup default observations of the component
        self.default_obs = {'pressure': 101325.0,
                            'CO2saturation': 0.0,
                            'mass_CO2_reservoir': 0.0}

        # Add accumulator observation mass_CO2_reservoir used within the component
        self.add_accumulator('mass_CO2_reservoir', sim=0.0)
        self.add_accumulator('volume_CO2_reservoir', sim=0.0)

        debug_msg = 'SimpleReservoir component created with name {}'.format(self.name)
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
            warn_msg = ''.join([
                'Parameter {} of SimpleReservoir component {} ',
                'is out of boundaries.']).format(key, self.name)
            if key.startswith('shale') and key.endswith('Thickness'):
                if (val < self.pars_bounds['shaleThickness'][0]) or (
                        val > self.pars_bounds['shaleThickness'][1]):
                    logging.warning(warn_msg)
            elif key.startswith('aquifer') and key.endswith('Thickness'):
                if (val < self.pars_bounds['aquiferThickness'][0]) or (
                        val > self.pars_bounds['aquiferThickness'][1]):
                    logging.warning(warn_msg)
            elif key in self.pars_bounds:
                if (val < self.pars_bounds[key][0]) or (
                        val > self.pars_bounds[key][1]):
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of SimpleReservoir component {}.']).format(key, self.name)
                logging.warning(warn_msg)


    def simulation_model(self, p, time_point=365.25, time_step=365.25,
                         injX=None, injY=None, locX=None, locY=None):
        """
        Return pressure and CO2 saturation at the bottom of leaking well.
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this modules docstring.

        :param p: input parameters of simple reservoir model
        :type p: dict

        :param time_point: time point (in days) for which the pressure and
            saturation are to be calculated; by default, its value is 365.25
            (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param mass_CO2_reservoir: mass of |CO2| accumulated in reservoir at
            time specified by time_point parameter; by default, its value is
            0.0; this parameter supplies value for an accumulator
            (needed for the next time point) if solution has to be calculated
            for several time steps.
        :type pressure: float

        :param injX: x-coordinate of the injector location
        :type injX: float

        :param injY: y-coordinate of the injector location
        :type injY: float

        :param locX: x-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locX: float

        :param locY: y-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locY: float

        :returns: out - dictionary of observations of simple reservoir
            component model;
            keys: ['pressure','CO2saturation','mass_CO2_reservoir']

        """
        debug_msg = ''.join([
            'SimpleReservoir component {} model call input: ',
            'parameters {}, location xy ({}, {})']).format(
                self.name, p, self.locX, self.locY)
        logging.debug(debug_msg)

        # Set coordinates of injector and leak point/point of interest
        if injX is None:
            injX = self.injX
        if injY is None:
            injY = self.injY
        if locX is None:
            locX = self.locX
        if locY is None:
            locY = self.locY

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Create instance of Parameters class
        inputParameters = sresrom.Parameters()

        # Set up number of shale layers
        nSL = int(actual_p['numberOfShaleLayers'])
        inputParameters.numberOfShaleLayers = nSL

        # Initialize shale and aquifer thickness arrays
        inputParameters.shaleThickness = actual_p['shaleThickness']*np.ones(nSL)
        inputParameters.aquiferThickness = actual_p['aquiferThickness']*np.ones((nSL-1))

        # Set up shale and aquifer thickness
        for i in range(nSL):
            nm = 'shale{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.shaleThickness[i] = p[nm]
        for i in range(nSL-1):
            nm = 'aquifer{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.aquiferThickness[i] = p[nm]

        inputParameters.reservoirThickness = actual_p['reservoirThickness']
        if inputParameters.reservoirThickness < 0:
            warn_msg = 'Reservoir thickness {0} less than expected lowest value of 0'.format(
                inputParameters.reservoirThickness)
            logging.warning(warn_msg)

        # Set up reservoir permeability
        inputParameters.reservoirPermeability = 10**actual_p['logResPerm']

        # Set up reservoir porosity
        inputParameters.reservoirPorosity = actual_p['reservoirPorosity']

        # Set up land surface pressure
        inputParameters.datumPressure = actual_p['datumPressure']

        # Set up brine and CO2 density
        inputParameters.brineDensity = actual_p['brineDensity']
        inputParameters.CO2Density = actual_p['CO2Density']

        # Set up brine and CO2 viscosity
        inputParameters.brineViscosity = actual_p['brineViscosity']
        inputParameters.CO2Viscosity = actual_p['CO2Viscosity']

        # Set up residual saturation and compressibility
        inputParameters.brineResidualSaturation = actual_p['brineResSaturation']
        inputParameters.compressibility = actual_p['compressibility']

        # Set up distance between injection and leaking well
        inputParameters.distance = sqrt((locX-injX)**2+(locY-injY)**2)

        # Rate at the injection well
        inputParameters.rate = actual_p['injRate']

        # Parameters from the class attributes
        inputParameters.timeStep = time_step       # in days
        inputParameters.timePoint = time_point     # in days
        inputParameters.prevCO2Volume = self.accumulators['volume_CO2_reservoir'].sim

        # Create solution object with defined input parameters
        sol = sresrom.Solution(inputParameters)

        # Initiate dictionary of leakage rates
        out = dict()
        if time_point == 0.0:   #  return initial state of the reservoir component
            sol.setup_initial_conditions()
            # Save observations
            out['pressure'] = sol.initialTopPressure
            out['CO2saturation'] = sol.interface/sol.thicknessH
        else:
            # Find solution corresponding to the inputParameters
            sol.find()
            # Save observations
            out['pressure'] = sol.pressure
            out['CO2saturation'] = sol.saturation

        out['mass_CO2_reservoir'] = sol.CO2Volume*inputParameters.CO2Density
        self.accumulators['mass_CO2_reservoir'].sim = sol.CO2Volume*inputParameters.CO2Density
        self.accumulators['volume_CO2_reservoir'].sim = sol.CO2Volume

        debug_msg = 'SimpleReservoir component {} model call output: {}'.format(
            self.name, out)
        logging.debug(debug_msg)

        # Return dictionary of outputs
        return out

    # Attributes for system connections
    needs_locXY = True
    needs_injXY = True
    system_params = ['numberOfShaleLayers',
                     'shale1Thickness',
                     'shale2Thickness',
                     'shale3Thickness',
                     'aquifer1Thickness',
                     'aquifer2Thickness',
                     'reservoirThickness',
                     'datumPressure']

    def reset(self):
        pass


def test_simple_reservoir_component():
    __spec__ = None

    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = iam_bc.SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('injRate', min=0.4, max=0.6, value=0.5)
    sres.add_par('shale1Thickness', min=30.0, max=50., value=40.0)
    sres.add_par('shale2Thickness', min=40.0, max=60., value=50.0)

    # Add observations of reservoir component model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with names sres.obsnm_0, sres.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Run system model using current values of its parameters
    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Print pressure and saturation
    # Since the observations at the particular time points are different variables,
    # method collect_observations_as_time_series creates lists of
    # values of observations belonging to a given component (e.g. cw) and having the same
    # common name (e.g. 'pressure', 'CO2saturation', etc) but differing in indices.
    # More details are given in the docstring and documentation to the method
    # collect_observations_as_time_series of SystemModel class.
    print('Pressure',
          sm.collect_observations_as_time_series(sres, 'pressure'), sep='\n')
    print('CO2 saturation',
          sm.collect_observations_as_time_series(sres, 'CO2saturation'), sep='\n')
    print('Mass of CO2 in reservoir',
          sm.collect_observations_as_time_series(sres, 'mass_CO2_reservoir'), sep='\n')

    # print('------------------------------------------------------------------')
    # print('                          UQ illustration ')
    # print('------------------------------------------------------------------')

    # import random
    # num_samples = 50
    # ncpus = 5
    # # Draw Latin hypercube samples of parameter values
    # s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))   # create sample set

    # # Run model using values in samples for parameter values
    # results = s.run(cpus=ncpus, verbose=False)

    # # Extract results from stochastic simulations
    # pressure = np.ones((num_samples, len(time_array)))
    # CO2_saturation = np.ones((num_samples, len(time_array)))
    # mass_CO2_reservoir = np.ones((num_samples, len(time_array)))

    # for ind in range(len(time_array)):
    #     pressure[:, ind] = s.recarray['sres.pressure_{}'.format(ind)]
    #     CO2_saturation[:, ind] = s.recarray['sres.CO2saturation_{}'.format(ind)]
    #     mass_CO2_reservoir[:, ind] = s.recarray['sres.mass_CO2_reservoir_{}'.format(ind)]

    # # Plot results
    # # Setup plot parameters
    # line_width = 1
    # xy_label_size = 14
    # title_size = 18
    # xlims = [0, 50]
    # fig = plt.figure(figsize=(18, 4))
    # ax = fig.add_subplot(131)
    # for j in range(num_samples):
    #     plt.plot(time_array/365.25, pressure[j]/1.0e+6,
    #              color="maroon", linewidth=line_width)
    # plt.xlabel('Time, t (years)', fontsize=xy_label_size)
    # plt.ylabel('Pressure, P (MPa)', fontsize=xy_label_size)
    # plt.title('Pressure: leaking well', fontsize=title_size)
    # plt.tight_layout()
    # plt.tick_params(labelsize=12)
    # plt.xlim(xlims)
    # plt.ylim([20, 45])
    # ax.get_yaxis().set_label_coords(-0.12, 0.5)

    # ax = fig.add_subplot(132)
    # for j in range(num_samples):
    #     plt.plot(time_array/365.25, CO2_saturation[j],
    #              color="green", linewidth=line_width)
    # plt.xlabel('Time, t (years)', fontsize=xy_label_size)
    # plt.ylabel('Saturation, S (-)', fontsize=xy_label_size)
    # plt.title(r'CO$_2$ saturation: leaking well', fontsize=title_size)
    # plt.tight_layout()
    # plt.tick_params(labelsize=12)
    # plt.xlim(xlims)
    # plt.ylim([0.0, 1.0])
    # ax.get_yaxis().set_label_coords(-0.12, 0.5)

    # ax = fig.add_subplot(133)
    # for j in range(num_samples):
    #     plt.plot(time_array/365.25, 479.0*mass_CO2_reservoir[j],
    #              color="darkblue", linewidth=line_width)
    # plt.xlabel('Time, t (years)', fontsize=xy_label_size)
    # plt.ylabel(r'Mass of CO$_2$, m (kg)', fontsize=xy_label_size)
    # plt.title(r'Amount of CO$_2$ in reservoir', fontsize=title_size)
    # plt.tight_layout()
    # plt.tick_params(labelsize=12)
    # plt.xlim(xlims)
    # plt.ylim([0.0, 5.0e+11])
    # ax.get_yaxis().set_label_coords(-0.12, 0.5)

    # to_save = False
    # if to_save:
    #     plt.savefig('ReservoirComponentCombinedPlot2.png', dpi=300)
