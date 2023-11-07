"""
Example illustrates linking of two analytical reservoir and multisegmented
wellbore components. The locations of the wellbores are random. Example also
shows Latin Hypercube sampling study.

Example of run:
$ python iam_sys_reservoir_mswell_2locs.py
"""

import sys
import os
import random

import matplotlib.pyplot as plt
import numpy as np

from openiam.components.iam_base_classes import SystemModel
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.multisegmented_wellbore_component import MultisegmentedWellbore
from openiam.matk import pyDOE


if __name__=='__main__':
    # For multiprocessing in Spyder
    __spec__ = None

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    xymins = np.array([100., 450.])
    xymaxs = np.array([300., 500.])
    num_wells = 2
    well_xys = xymins + pyDOE.lhs(2, samples=num_wells)*(xymaxs-xymins)

    aress = []
    mss = []
    for i, crds in enumerate(well_xys):

        # Add reservoir components
        aress.append(sm.add_component_model_object(
            AnalyticalReservoir(name='ares'+str(i), parent=sm,
                                injX=0., injY=0., locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        aress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        aress[-1].add_par('injRate', value=0.8, vary=False)
        aress[-1].add_par('shale1Thickness', min=30.0, max=50., value=40.0)
        aress[-1].add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
        aress[-1].add_par('aquifer2Thickness', min=30.0, max=50., value=45.0)
        aress[-1].add_par('logResPerm', min=-13.5, max=-12.5, value=-13.)
        aress[-1].add_par('reservoirRadius', value=1000., vary=False)

        # Add observations of reservoir component model to be used by the next component
        aress[-1].add_obs_to_be_linked('pressure')
        aress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        aress[-1].add_obs('pressure')
        aress[-1].add_obs('CO2saturation')

        # Add multisegmented wellbore components
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name='ms'+str(i), parent=sm)))

        # A lot of parameters of multisegmented wellbore component
        # are the same as for the reservoir component
        # Add parameters linked to the same parameters from reservoir model
        mss[-1].add_par_linked_to_par(
            'numberOfShaleLayers', aress[-1].deterministic_pars['numberOfShaleLayers'])
        mss[-1].add_par_linked_to_par(
            'shale1Thickness', aress[-1].pars['shale1Thickness'])
        mss[-1].add_par_linked_to_par(
            'shale2Thickness', aress[-1].default_pars['shaleThickness'])
        mss[-1].add_par_linked_to_par(
            'shale3Thickness', aress[-1].default_pars['shaleThickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer1Thickness', aress[-1].pars['aquifer1Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer2Thickness', aress[-1].pars['aquifer2Thickness'])
        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_kwarg_linked_to_obs('pressure', aress[-1].linkobs['pressure'])
        mss[-1].add_kwarg_linked_to_obs('CO2saturation', aress[-1].linkobs['CO2saturation'])
        mss[-1].add_obs('brine_aquifer1')
        mss[-1].add_obs('CO2_aquifer1')

    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    linespec = ['-r', '-b', '-g', '-k', '-m']
    f1, ax = plt.subplots(2, 2, figsize=(20, 12))
    # Print and plot pressure and saturation
    for i, ares in enumerate(aress):
        pressure = sm.collect_observations_as_time_series(ares, 'pressure')
        CO2saturation = sm.collect_observations_as_time_series(ares, 'CO2saturation')
        print('Pressure at wellbore {}:'.format(i+1), pressure, sep='\n')
        print('CO2 saturation at wellbore {}:'.format(i+1), CO2saturation, sep='\n')
        print('------------------------------------------------------------------')
        ax[0, 0].plot(time_array/365.25, pressure, linespec[i], label='wellbore '+str(i+1))
        ax[0, 0].set_xlabel('Time, [years]')
        ax[0, 0].set_ylabel('Pressure, [Pa]')
        ax[0, 1].plot(time_array/365.25, CO2saturation, linespec[i], label='wellbore '+str(i+1))
        ax[0, 1].set_xlabel('Time, [years]')
        ax[0, 1].set_ylabel(r'CO$_2$ saturation, [-]')

    for i, ms in enumerate(mss):
        CO2aq1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
        brineaq1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')
        print('Leakage rates at wellbore {}:'.format(i+1))
        print('CO2: aquifer 1', CO2aq1, sep='\n')
        print('Brine: aquifer 1', brineaq1, sep='\n')
        ax[1, 0].plot(time_array/365.25, CO2aq1, linespec[i], label='wellbore '+str(i+1))
        ax[1, 0].set_xlabel('Time, [years]')
        ax[1, 0].set_ylabel(r'CO$_2$ leakage rates to aquifer 1, [kg/s]')
        ax[1, 1].plot(time_array/365.25, brineaq1, linespec[i], label='wellbore '+str(i+1))
        ax[1, 1].set_xlabel('Time, [years]')
        ax[1, 1].set_ylabel('Brine leakage rates to aquifer 1, [kg/s]')

    # Setup legend on sublots
    ax[0, 0].legend()
    ax[0, 1].legend()
    ax[1, 0].legend()
    ax[1, 1].legend()
    plt.show()

    print('------------------------------------------------------------------')
    print('                          UQ illustration ')
    print('------------------------------------------------------------------')
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=30, seed=random.randint(500, 1100))   # create sample set

    # Run model using values in samples for parameter values
    results = s.run(cpus=4, verbose=False)

    print('Results of simulations')
    print(results)
