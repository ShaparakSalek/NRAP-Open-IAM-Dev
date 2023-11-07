'''
This example couples the analytical reservoir, multisegmented wellbore and
carbonate aquifer models. The saturation/pressure output produced by analytical
reservoir model is used to drive leakage from a single multisegmented wellbore
model, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates and cumulative mass fluxes to the
carbonate aquifer models representing 2 aquifers in the system.

Example of run:
$ python iam_sys_reservoir_mswell_2aquifers.py
'''

import sys
import os
import numpy as np

from openiam.components.iam_base_classes import SystemModel
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.multisegmented_wellbore_component import MultisegmentedWellbore
from openiam.components.rate_to_mass_adapter import RateToMassAdapter
from openiam.components.carbonate_aquifer_component import CarbonateAquifer


if __name__ == "__main__":
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('numberOfShaleLayers', value=3, vary=False)
    ares.add_par('shale1Thickness', min=10.0, max=25., value=15.0)
    ares.add_par('shale2Thickness', min=10.0, max=25., value=25.0)
    ares.add_par('shale3Thickness', min=300.0, max=500., value=350.0)
    ares.add_par('injRate', min=0.1, max=10.0, value=1.0)
    ares.add_par('aquifer1Thickness', min=100.0, max=150., value=120.0)
    ares.add_par('aquifer2Thickness', min=100.0, max=140., value=110.0)
    ares.add_par('reservoirThickness', min=30.0, max=50., value=40.0)
    ares.add_par('logResPerm', min=-13.5, max=-12.5, value=-13)
    ares.add_par('injRate', value=0.025, vary=False)

    # Add observations of reservoir component model
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(
        MultisegmentedWellbore(name='ms', parent=sm))
    ms.add_par('wellRadius', min=0.01, max=0.02, value=0.015)
    ms.add_par('logWellPerm', min=-14.0, max=-12.0, value=-12.5)

    # Add linked parameters: common to reservoir and wellbore components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             ares.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness', ares.pars['shale1Thickness'])
    ms.add_par_linked_to_par('shale2Thickness', ares.pars['shale2Thickness'])
    ms.add_par_linked_to_par('shale3Thickness', ares.pars['shale3Thickness'])
    ms.add_par_linked_to_par('aquifer1Thickness', ares.pars['aquifer1Thickness'])
    ms.add_par_linked_to_par('aquifer2Thickness', ares.pars['aquifer2Thickness'])
    ms.add_par_linked_to_par('reservoirThickness', ares.pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure',
                             ares.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs_to_be_linked('CO2_aquifer1')
    ms.add_obs_to_be_linked('CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_aquifer1')
    ms.add_obs_to_be_linked('brine_aquifer2')
    ms.add_obs_to_be_linked('mass_CO2_aquifer1')
    ms.add_obs_to_be_linked('mass_CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs_to_be_linked('CO2_atm')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')

    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(
        RateToMassAdapter(name='adapt', parent=sm))
    adapt.add_kwarg_linked_to_collection('CO2_aquifer1',
        [ms.linkobs['CO2_aquifer1'], ms.linkobs['CO2_aquifer2']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer2',
        [ms.linkobs['CO2_aquifer2'], ms.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer1',
        [ms.linkobs['brine_aquifer1'], ms.linkobs['brine_aquifer2']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer2',
        [ms.linkobs['brine_aquifer2'], ms.linkobs['brine_atm']])
    adapt.add_obs_to_be_linked('mass_CO2_aquifer1')
    adapt.add_obs_to_be_linked('mass_CO2_aquifer2')
    adapt.add_obs_to_be_linked('mass_brine_aquifer1')
    adapt.add_obs_to_be_linked('mass_brine_aquifer2')
    adapt.add_obs('mass_CO2_aquifer1')
    adapt.add_obs('mass_brine_aquifer1')
    adapt.add_obs('mass_CO2_aquifer2')
    adapt.add_obs('mass_brine_aquifer2')

    # Add 2 carbonate aquifer model objects and define parameters
    cas = []
    for ind in range(2):
        cas.append(sm.add_component_model_object(
            CarbonateAquifer(name='ca'+str(ind),parent=sm)))
        cas[-1].add_par('perm_var', min=0.017, max=1.89, value=0.9535)
        cas[-1].add_par('corr_len', min=1.0, max=3.95, value=2.475)
        cas[-1].add_par('aniso', min=1.1, max=49.1, value=25.1)
        cas[-1].add_par('mean_perm', min=-13.8, max=-10.3, value=-12.05)
        cas[-1].add_par_linked_to_par(
            'aqu_thick', ares.pars['aquifer{}Thickness'.format(ind+1)])
        cas[-1].add_par('hyd_grad', min=2.88e-4, max=1.89e-2, value=9.59e-3)
        cas[-1].add_par('calcite_ssa', min=0, max=1.e-2, value=5.5e-03)
        cas[-1].add_par('organic_carbon', min=0, max=1.e-2, value=5.5e-03)
        cas[-1].add_par('benzene_kd', min=1.49, max=1.73, value=1.61)
        cas[-1].add_par('benzene_decay', min=0.15, max=2.84, value=1.5)
        cas[-1].add_par('nap_kd', min=2.78, max=3.18, value=2.98)
        cas[-1].add_par('nap_decay', min=-0.85, max=2.04, value=0.595)
        cas[-1].add_par('phenol_kd', min=1.21, max=1.48, value=1.35)
        cas[-1].add_par('phenol_decay', min=-1.22, max=2.06, value=0.42)
        cas[-1].add_par('cl', min=0.1, max=6.025, value=0.776)
        cas[-1].model_kwargs['x'] = [0.]
        cas[-1].model_kwargs['y'] = [0.]

        CO2_rate_obs_list = []
        brine_rate_obs_list = []
        CO2_mass_obs_list = []
        brine_mass_obs_list = []
        CO2_rate_obs_list.append(ms.linkobs['CO2_aquifer'+str(ind+1)])
        brine_rate_obs_list.append(ms.linkobs['brine_aquifer'+str(ind+1)])
        CO2_mass_obs_list.append(adapt.linkobs['mass_CO2_aquifer'+str(ind+1)])
        brine_mass_obs_list.append(adapt.linkobs['mass_brine_aquifer'+str(ind+1)])

        # Add aquifer component's keyword argument co2_rate linked to the collection created above
        cas[-1].add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
        cas[-1].add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)
        cas[-1].add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
        cas[-1].add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

        # Add observations (output) from the carbonate aquifer model
        cas[-1].add_obs('pH_volume')
        cas[-1].add_obs('TDS_volume')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('pressure',
          sm.collect_observations_as_time_series(ares, 'pressure'))
    print('------------------------------------------------------------------')
    print('CO2saturation',
          sm.collect_observations_as_time_series(ares, 'CO2saturation'))
    print('------------------------------------------------------------------')
    print('CO2_aquifer1',
          sm.collect_observations_as_time_series(ms, 'CO2_aquifer1'))
    print('------------------------------------------------------------------')
    print('brine_aquifer1',
          sm.collect_observations_as_time_series(ms, 'brine_aquifer1'))
    print('------------------------------------------------------------------')
    print('mass_CO2_aquifer1',
          sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer1'))
    print('------------------------------------------------------------------')
    print('mass_brine_aquifer1',
          sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer1'))
    print('------------------------------------------------------------------')
    print('mass_CO2_aquifer2',
          sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer2'))
    print('------------------------------------------------------------------')
    print('mass_brine_aquifer2',
          sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer2'))
    print('------------------------------------------------------------------')
    print('aquifer1: pH',
          sm.collect_observations_as_time_series(cas[0], 'pH_volume'))
    print('------------------------------------------------------------------')
    print('aquifer1: TDS',
          sm.collect_observations_as_time_series(cas[0], 'TDS_volume'))
    print('------------------------------------------------------------------')
    print('aquifer2: pH',
          sm.collect_observations_as_time_series(cas[1], 'pH_volume'))
    print('------------------------------------------------------------------')
    print('aquifer2: TDS',
          sm.collect_observations_as_time_series(cas[1], 'TDS_volume'))
