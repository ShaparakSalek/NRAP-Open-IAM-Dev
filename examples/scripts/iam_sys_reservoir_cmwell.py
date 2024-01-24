'''
System model contains two linked component models: reservoir model and
cemented wellbore model. Some of the parameters are deterministic, some are
stochastic, some are composite or linked to other parameters.
Keyword arguments of the second component are obtained from the observations
of the first component. Plots of leakage rates are returned.

Examples of run:
$ python iam_sys_reservoir_cmwell.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

from openiam.components.iam_base_classes import SystemModel
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.cemented_wellbore_component import CementedWellbore


if __name__=='__main__':
    # Create system model
    num_years = 50.
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('numberOfShaleLayers', value=3, vary=False)
    ares.add_par('shale1Thickness', min=500.0, max=550., value=525.0)
    ares.add_par('shale2Thickness', min=450.0, max=500., value=475.0)
    # Shale 3 has a fixed thickness of 11.2 m
    ares.add_par('shale3Thickness', value=11.2, vary=False)
    # Aquifer 1 (thief zone has a fixed thickness of 22.4)
    ares.add_par('aquifer1Thickness', value=22.4, vary=False)
    # Aquifer 2 (shallow aquifer) has a fixed thickness of 19.2
    ares.add_par('aquifer2Thickness', value=19.2, vary=False)
    # Reservoir has a fixed thickness of 51.2
    ares.add_par('reservoirThickness', value=51.2, vary=False)
    ares.add_par('injRate', value=0.025, vary=False)

    # Add observations of reservoir component model
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')

    # Add cemented wellbore component
    cw = sm.add_component_model_object(CementedWellbore(name='cw', parent=sm))

    # Add parameters of cemented wellbore component
    cw.add_par('logWellPerm', min=-13.9, max=-12.0, value=-12.0)

    # Add keyword arguments of the cemented wellbore component model
    cw.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    cw.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])

    # Add composite parameters of cemented wellbore component
    # Here, we illustrate two ways to define expressions for composite parameters
    # One way
    cw.add_composite_par('wellDepth', expr=ares.pars['shale1Thickness'].name+
        '+'+ares.pars['shale2Thickness'].name+
        '+'+ares.deterministic_pars['shale3Thickness'].name+
        '+'+ares.deterministic_pars['aquifer1Thickness'].name+
        '+'+ares.deterministic_pars['aquifer2Thickness'].name)
    # Second shorter (and more explicit) way
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

    # Run forward simulation
    sm.forward()

    # Collect observations
    CO2_leakrates_aq1 = sm.collect_observations_as_time_series(cw, 'CO2_aquifer1')
    CO2_leakrates_aq2 = sm.collect_observations_as_time_series(cw, 'CO2_aquifer2')

    brine_leakrates_aq1 = sm.collect_observations_as_time_series(cw, 'brine_aquifer1')
    brine_leakrates_aq2 = sm.collect_observations_as_time_series(cw, 'brine_aquifer2')

    # Print results: CO2 and brine leakage rates and pressure/saturation at the wellbore
    print('CO2 leakage rates to aquifer 1:', CO2_leakrates_aq1, sep='\n')
    print('CO2 leakage rates to aquifer 2:', CO2_leakrates_aq2, sep='\n')
    print('Brine leakage rates to aquifer 1:', brine_leakrates_aq1, sep='\n')
    print('Brine leakage rates to aquifer 2:', brine_leakrates_aq2, sep='\n')
    print('Pressure:', sm.collect_observations_as_time_series(ares,'pressure'), sep='\n')
    print('CO2 saturation:', sm.collect_observations_as_time_series(ares, 'CO2saturation'), sep='\n')

    # Plot CO2 and brine leakage rates along the wellbore
    plt.figure(1)
    plt.plot(sm.time_array/365.25, CO2_leakrates_aq1, color='#000055',
             linewidth=2, label="aquifer 1")
    plt.plot(sm.time_array/365.25, CO2_leakrates_aq2, color='#FF0066',
             linewidth=2, label="aquifer 2")
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.xlabel('Time, t [years]')
    plt.ylabel('Leakage rates, q [kg/s]')
    plt.title(r'Leakage of CO$_2$ to aquifer 1 and aquifer 2')

    plt.figure(2)
    plt.plot(sm.time_array/365.25, brine_leakrates_aq1, color='#000055',
             linewidth=2, label="aquifer 1")
    plt.plot(sm.time_array/365.25, brine_leakrates_aq2, color='#FF0066',
             linewidth=2, label="aquifer 2")
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.xlabel('Time, t [years]')
    plt.ylabel('Leakage rates, q [kg/s]')
    plt.title(r'Leakage of brine to aquifer 1 and aquifer 2')
    plt.show()
