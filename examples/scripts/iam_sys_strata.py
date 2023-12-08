"""
Example demonstrating the use of the Stratigraphy component. Here, the component 
is given unit thicknesses and set up with all of the composite depth parameters 
for a Stratigraphy component. These composite parameters represent the bottom, 
middle, or top of each unit in the domain (shales, aquifers, and the reservoir). 
These names and expressions (strings used to calcualte the values) for these 
composite parameters can be set using the get_composite_depth_names() and 
get_depth_expr() methods of the Stratigraphy component.

After the simulation has been run, the value of each composite depth parameter 
is printed.

Example of run:
$ python iam_sys_strata.py
"""

import sys
import os
import logging
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, Stratigraphy


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    # Define keyword arguments of the system model. The times used do not have 
    # an impact on the composite depth parameters printed further below.
    num_years = 10
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    
    numberOfShaleLayers = 5
    reservoir_thickness = 31
    shale_thicknesses = [322, 216, 67, 156, 34]
    aquifer_thicknesses = [51, 45, 142, 12]
    
    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)
    
    # Add stratigraphy component
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))
    
    strata.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    strata.add_par('reservoirThickness', value=reservoir_thickness, vary=False)
    
    for unit_num in range(1, numberOfShaleLayers + 1):
        strata.add_par('shale{}Thickness'.format(unit_num), 
                       value=shale_thicknesses[unit_num - 1], vary=False)
        
        if unit_num < numberOfShaleLayers:
            strata.add_par('aquifer{}Thickness'.format(unit_num), 
                           value=aquifer_thicknesses[unit_num - 1], vary=False)
    
    composite_depth_pars = strata.get_composite_depth_names()
    
    for comp_depth in composite_depth_pars:
        par_expr = strata.get_depth_expr(comp_depth)
        strata.add_composite_par(comp_depth, par_expr)
    
    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    
    # Collect composite depth parameter values. The composite parameters are 
    # calculated during the simulation. If you tried to print these values 
    # without first running the simulation, all of the values would be zero.
    print('Depths (m) from composite parameters:')
    for comp_depth in sorted(composite_depth_pars):
        comp_par_val = sm.component_models['strata'].composite_pars[
            comp_depth].value
        print(comp_depth + ': ', comp_par_val)
