'''
Example illustrates system model containing only analytical reservoir model.
System model is run for an array of time points. Parameters of the reservoir
model are defined as default, deterministic and stochastic. Plots of pressure
and |CO2| saturation are returned as output.

Examples of run:
$ python iam_sys_reservoir.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

from openiam.components.iam_base_classes import SystemModel
from openiam.components.analytical_reservoir_component import AnalyticalReservoir


if __name__=='__main__':

    # Create system model
    time_array = np.arange(0, 201)   # in days
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('numberOfShaleLayers', value=3, vary=False)
    ares.add_par('injRate', min=0.1, max=0.4, value=0.2)
    ares.add_par('shale1Thickness', min=25.0, max=80., value=50.0)
    ares.add_par('shale2Thickness', min=35.0, max=90., value=45.0)

    # Add observations of reservoir component model
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')

    # Run system model using current values of its parameters
    sm.forward()

    # Assign observations of the model to pressure and CO2saturation variables
    pressure = sm.collect_observations_as_time_series(ares, 'pressure')
    CO2saturation = sm.collect_observations_as_time_series(ares, 'CO2saturation')

    # Print pressure and saturation
    print('pressure', pressure, sep='\n')
    print('CO2 saturation', CO2saturation, sep='\n')

    # Plot pressure and saturation
    plt.figure(1)
    plt.plot(time_array, pressure/1.0e+6,
             'c-', linewidth=2, label='pressure')
    plt.xlabel('Time, t (years)')
    plt.ylabel('Pressure, P (MPa)')

    plt.figure(2)
    plt.plot(time_array, CO2saturation,
             'r-', linewidth=2, label='saturation')
    plt.xlabel('Time, t (years)')
    plt.ylabel('Saturation, S')
