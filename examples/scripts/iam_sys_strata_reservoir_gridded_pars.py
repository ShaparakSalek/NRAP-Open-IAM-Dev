# -*- coding: utf-8 -*-
"""
This example couples stratigraphy and analytical reservoir components.
Stratigraphy component defines gridded parameter shale1Thickness which is
then used to define a deterministic parameter of the same name
for the analytical reservoir component.

Example of run:
$ python iam_sys_strata_reservoir_gridded_pars.py
"""
import os
import sys
import logging
import numpy as np

from openiam.components.iam_base_classes import SystemModel
from openiam.components.stratigraphy_component import Stratigraphy
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.iam_gridded_observation import DataInterpolator


if __name__ == "__main__":
    # Define logging level
    logging.basicConfig(level=logging.WARNING)

    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.arange(0.0,num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    int1 = DataInterpolator(name='int1', parent=sm,
                            header_file_dir=os.path.join(
                                '..', '..', 'data', 'reservoir',
                                'lookuptables', 'Test_2d'),
                            data_file='shale1Thickness.csv')
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of reservoir component model
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_gridded_par('shale1Thickness', int1)
    strata.add_par('shale2Thickness', min=40.0, max=60., value=50.0)

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm,
        injX=37200., injY=48000., locX=37478.0, locY=48333.0))  # chosen injX, injY are arbitrary

    # Add parameters of reservoir component model
    ares.add_par_linked_to_par('numberOfShaleLayers',
                         strata.deterministic_pars['numberOfShaleLayers'])
    ares.add_par('injRate', min=0.4, max=0.6, value=0.5)
    ares.add_par('reservoirRadius', value=500.0, vary=False)
    ares.add_par('shale1Thickness', vary=False,
        value=strata.gridded_pars['shale1Thickness'](np.array([[ares.locX, ares.locY]])))
    ares.add_par_linked_to_par('shale2Thickness', strata.pars['shale2Thickness'])
    ares.add_par_linked_to_par('shale3Thickness', strata.default_pars['shaleThickness'])

    ares.add_par_linked_to_par('aquifer1Thickness', strata.default_pars['aquiferThickness'])
    ares.add_par_linked_to_par('aquifer2Thickness', strata.default_pars['aquiferThickness'])
    ares.add_par_linked_to_par('aquifer3Thickness', strata.default_pars['aquiferThickness'])

    ares.add_par_linked_to_par('reservoirThickness', strata.default_pars['reservoirThickness'])

    ares.add_par_linked_to_par('datumPressure', strata.default_pars['datumPressure'])

    ares.add_obs('pressure')

    # Run system model using current values of its parameters
    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Print pressure
    print('pressure', sm.collect_observations_as_time_series(ares, 'pressure'), sep='\n')
