# -*- coding: utf-8 -*-
'''
This example couples the analytical reservoir, open wellbore and
atmosphere models. The saturation/pressure output produced by analytical
reservoir model is used to drive leakage from a single open wellbore
model, |CO2| leakage rates are passed to the atmosphere model.
Plots of relevant observations are created.

Example of run:
$ python iam_sys_reservoir_openwell_atmosphere.py
'''

import sys
import os

import numpy as np
import matplotlib.pyplot as plt

from openiam.components.iam_base_classes import SystemModel
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.open_wellbore_component import OpenWellbore
from openiam.components.rate_to_mass_adapter import RateToMassAdapter
from openiam.components.atmRom_component import AtmosphericROM


if __name__ == "__main__":

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}  # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('numberOfShaleLayers', value=3, vary=False)
    ares.add_par('shale1Thickness', min=900.0, max=1100., value=1000.0)
    ares.add_par('shale2Thickness', min=900.0, max=1100., value=1000.0)
    # Shale 3 has a fixed thickness of 11.2 m
    ares.add_par('shale3Thickness', value=11.2, vary=False)
    # Aquifer 1 (thief zone has a fixed thickness of 22.4)
    ares.add_par('aquifer1Thickness', value=22.4, vary=False)
    # Aquifer 2 (shallow aquifer) has a fixed thickness of 19.2
    ares.add_par('aquifer2Thickness', value=500, vary=False)
    # Reservoir has a fixed thickness of 51.2
    ares.add_par('reservoirThickness', value=51.2, vary=False)

    # Add observations of reservoir component model
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')

    # Add open wellbore component
    ow = sm.add_component_model_object(OpenWellbore(name='ow', parent=sm))

    # Add parameters of open wellbore component
    ow.add_par('wellRadius', min=0.01, max=0.02, value=0.015)
    ow.add_par('logReservoirTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par('logAquiferTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par('brineSalinity', value=0.1, vary=False)
    ow.add_par('wellTop', value=0.0, vary=False)

    # Add keyword arguments of the open wellbore component model
    ow.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    ow.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])

    # Add composite parameter of open wellbore component
    ow.add_composite_par('reservoirDepth',
                         expr='+'.join(['ares.shale1Thickness',
                                        'ares.shale2Thickness',
                                        'ares.shale3Thickness',
                                        'ares.aquifer1Thickness',
                                        'ares.aquifer2Thickness']))

    # Add observations of open wellbore component model
    ow.add_obs_to_be_linked('CO2_atm')
    ow.add_obs('CO2_aquifer')  # zero since well top is in aquifer
    ow.add_obs('brine_aquifer')  # zero since well top is in aquifer
    ow.add_obs('CO2_atm')
    ow.add_obs('brine_atm')

    # Add Atm ROM component
    satm = sm.add_component_model_object(AtmosphericROM(name='satm', parent=sm))
    satm.add_par('T_amb', value=20.0)
    satm.add_par('P_amb', value=1.01E+00)
    satm.add_par('V_wind', value=5.0)
    satm.add_par('C0_critical', value=0.01)

    satm.model_kwargs['x_coor'] = [ares.locX]
    satm.model_kwargs['y_coor'] = [ares.locY]

    satm.model_kwargs['x_receptor'] = [100, 110, 200, 250]
    satm.model_kwargs['y_receptor'] = [110, 100, 100, 100]

    # Create collector for leakrates
    co2_leakrates_collector = []
    co2_leakrates_collector.append(ow.linkobs['CO2_atm'])

    satm.add_kwarg_linked_to_collection('co2_leakrate', co2_leakrates_collector)

    # Add observations for receptors
    for i, r in enumerate(satm.model_kwargs['x_receptor']):
        satm.add_obs('outflag_r{0:03}'.format(i))

    n_sources = len(satm.model_kwargs['x_coor'])
    for i in range(n_sources):
        satm.add_obs('x_new_s{0:03}'.format(i))
        satm.add_obs('y_new_s{0:03}'.format(i))
        satm.add_obs('critical_distance_s{0:03}'.format(i))

    satm.add_obs('num_sources')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('pressure',
          sm.collect_observations_as_time_series(ares, 'pressure'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2saturation',
          sm.collect_observations_as_time_series(ares, 'CO2saturation'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2_aquifer',
          sm.collect_observations_as_time_series(ow, 'CO2_aquifer'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2_atm',
          sm.collect_observations_as_time_series(ow, 'CO2_atm'), sep='\n')
    print('------------------------------------------------------------------')
    num_sources = sm.collect_observations_as_time_series(satm, 'num_sources')
    print('num_sources', num_sources)
    for i, r in enumerate(satm.model_kwargs['x_receptor']):
        print('------------------------------------------------------------------')
        print('Receptor_{}'.format(i),
              sm.collect_observations_as_time_series(
                  satm, 'outflag_r{0:03}'.format(i)), sep='\n')

    xnew = []
    ynew = []
    crit_dist = []
    for i in range(n_sources):
        xnew.append(sm.collect_observations_as_time_series(
            satm, 'x_new_s{0:03}'.format(i)))
        ynew.append(sm.collect_observations_as_time_series(
            satm, 'y_new_s{0:03}'.format(i)))
        crit_dist.append(sm.collect_observations_as_time_series(
            satm, 'critical_distance_s{0:03}'.format(i)))

    # For first source only
    plt.plot(time_array/365.25, crit_dist[0])
    plt.xlabel('Time $(y)$')
    plt.ylabel('Critical Distance from 1st source $(m)$')

    plt.show()

    print('------------------------------------------------------------------')
    # Limit to first 5 time steps for brevity
    message = 'At time {time} there is {num_sources} source(s):'
    for i, ns in enumerate(num_sources[:5]):
        print(message.format(time=time_array[i]/365.25, num_sources=ns))
        for j in range(ns):
            print('Point: ({x}, {y})  Critical_distance: {cd}'.format(
                x=xnew[j][i], y=ynew[j][i], cd=crit_dist[j][i]), sep='\n')
