'''
This example demonstrates the use of the DippingStratigraphy component and 
the Stratigraphy plot type. The stratigraphy_plot function in stratigraphy_plot.py 
requires the input dictionaries yaml_data and model_data. These dictionaries can 
be made with .yaml Control Files, but they can also be made manually as 
demonstrated here. Note that the function requires 
yaml_data['Plots'][plotName]['Stratigraphy'], where plotName is an arbitrary
user defined string (e.g., 'Plot1' or 'StratPlot'). The function does not
require any of the plotting options that can be contained within
yaml_data['Plots'][plotName]['Stratigraphy']. If any options are not present,
the default settings will be used.

Example of run:
$ python iam_sys_reservoir_mswell_stratplot_dipping_strata.py
'''

import sys
import os
import numpy as np
import random
import datetime
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, DippingStratigraphy, AnalyticalReservoir,
                     MultisegmentedWellbore)

import openiam as iam
import openiam.visualize as iam_vis


if __name__ == "__main__":
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Stratigraphy information
    datumPressure = 101325
    numberOfShaleLayers = 3
    shale1Thickness = 650
    shale2Thickness = 450
    shale3Thickness = 150
    aquifer1Thickness = 200
    aquifer2Thickness = 100
    reservoirThickness = 75

    # These lists are required by stratigraphy_plot()
    shaleThicknessesReferencePoint = [shale1Thickness, shale2Thickness,
                                      shale3Thickness]
    aquiferThicknessesReferencePoint = [aquifer1Thickness, aquifer2Thickness]

    strike = 315
    dip = 5
    dipDirection = 'NE'
    locXRef = 0
    locYRef = 0

    injectionX = 2500
    injectionY = 2500

    max_x_value = 5000
    max_y_value = 5000

    # Set up the dictionaries required by stratigraphy_plot()
    model_data = dict()
    model_data['OutputDirectory'] = os.path.join(
        iam.IAM_DIR, 'output', 'stratplot_example_dipping_strata_'
        + str(datetime.date.today()))

    yaml_data = dict()
    
    # A 'DippingStratigraphy' key needs to be in yaml_data - the TTFD plot type 
    # checks yaml_data to see what kind of stratigraphy is being used.
    yaml_data['DippingStratigraphy'] = dict()
    
    # yaml_data also needs a dictionary containing all of the parameters
    pars_for_yaml = {
        'numberOfShaleLayers': numberOfShaleLayers, 'datumPressure': datumPressure, 
        'shale1Thickness': shale1Thickness, 'shale2Thickness': shale2Thickness, 
        'shale3Thickness': shale3Thickness, 'aquifer1Thickness': aquifer1Thickness, 
        'aquifer2Thickness': aquifer2Thickness, 'reservoirThickness': reservoirThickness}
    
    yaml_data['DippingStratigraphy']['Parameters'] = pars_for_yaml
    
    plotName = 'Plot1.png'
    yaml_data['Plots'] = dict()
    yaml_data['Plots'][plotName] = dict()
    yaml_data['Plots'][plotName]['Stratigraphy'] = dict()

    use_specific_settings = True
    # I do this to demonstrate that the remaining options are not required. Try
    # setting use_specific_settings to False.
    if use_specific_settings:
        for key in ['PlotWellbores', 'PlotWellLabels', 'PlotStratComponents',
                    'PlotInjectionSites', 'PlotInjectionSiteLabels']:
            yaml_data['Plots'][plotName]['Stratigraphy'][key] = True

        yaml_data['Plots'][plotName]['Stratigraphy']['FigureDPI'] = 100

        yaml_data['Plots'][plotName]['Stratigraphy']['SaveCSVFiles'] = False

        yaml_data['Plots'][plotName]['Stratigraphy']['SpecifyXandYLims'] = \
            {'xLims': [-100, max_x_value + 100],
             'yLims': [-100, max_y_value + 100]}

        yaml_data['Plots'][plotName]['Stratigraphy']['SpecifyXandYGridLims'] = \
            {'gridXLims': [0, max_x_value], 'gridYLims': [0, max_y_value]}

        yaml_data['Plots'][plotName]['Stratigraphy']['View'] = \
            {'ViewAngleElevation': [10, 10, 10],
             'ViewAngleAzimuth': [300, 315, 330]}

        yaml_data['Plots'][plotName]['Stratigraphy']['StrikeAndDipSymbol'] = \
            {'PlotSymbol': True, 'coordx': 500, 'coordy': 1000, 'length': 500}

    # Create system model
    sm = SystemModel(model_kwargs = sm_model_kwargs)

    # Lists of components
    strata = []
    ares = []
    ms = []

    # List of all components except for stratigraphy components. This list is
    # needed for the stratigraphy_plot() function.
    components = []
    
    numberOfWells = 6
    np.random.seed(random.randint(0, 1.0e6))
    well_x_values = np.random.rand(1, numberOfWells) * max_x_value
    well_x_values = well_x_values.tolist()[0]
    well_y_values = np.random.rand(1, numberOfWells) * max_y_value
    well_y_values = well_y_values.tolist()[0]

    for locRef in range(0, numberOfWells):
        # The names for the reservoir and wellbore components at this location.
        aresName = 'AnalyticalReservoir_{0:03}'.format(locRef)
        msName = 'MultisegmentedWellbore_{0:03}'.format(locRef)
        # When spatially variable stratigraphy is used in the control file
        # interface, the stratigraphy component (here, DippingStratigraphy) made
        # for each component is named 'strata' + the component's name. Some
        # visualization codes, like stratigraphy_plot.py, are made to use this
        # naming convention.
        strataName = 'strata' + msName
        
        strata.append(sm.add_component_model_object(DippingStratigraphy(
            name=strataName, parent=sm, locXRef=locXRef, locYRef=locYRef, 
            dipDirection=dipDirection, locX=well_x_values[locRef], locY=well_y_values[locRef])))
        
        # Add parameters for the dipping stratrigraphy component
        strata[-1].add_par('numberOfShaleLayers',
                           value=numberOfShaleLayers, vary=False)
        strata[-1].add_par('shale1Thickness',
                           value=shale1Thickness, vary=False)
        strata[-1].add_par('shale2Thickness',
                           value=shale2Thickness, vary=False)
        strata[-1].add_par('shale3Thickness',
                           value=shale3Thickness, vary=False)
        strata[-1].add_par('aquifer1Thickness',
                           value=aquifer1Thickness, vary=False)
        strata[-1].add_par('aquifer2Thickness',
                           value=aquifer2Thickness, vary=False)
        strata[-1].add_par('reservoirThickness',
                           value=reservoirThickness, vary=False)
        strata[-1].add_par('datumPressure', value=datumPressure, vary=False)

        strata[-1].add_par('strike', value=strike, vary=False)
        strata[-1].add_par('dip', value=dip, vary=False)
        
        # Only use get_thickness_obs_names() and get_depth_obs_names() after the 
        # numberOfShaleLayers parameter has been assigned.
        thickness_obs = strata[-1].get_thickness_obs_names()
        
        # The thicness and depths observations are only produced in the first time 
        # step, so use index=[0] when adding these observations.
        for ob_nm in thickness_obs:
            strata[-1].add_obs(ob_nm, index=[0])
            strata[-1].add_obs_to_be_linked(ob_nm)
        
        depth_obs = strata[-1].get_depth_obs_names()
        
        for ob_nm in depth_obs:
            strata[-1].add_obs(ob_nm, index=[0])
            strata[-1].add_obs_to_be_linked(ob_nm)
        
        # Add reservoir component for the current well
        ares.append(sm.add_component_model_object(AnalyticalReservoir(
            name=aresName, parent=sm, locX=well_x_values[locRef], 
            locY=well_y_values[locRef], injX=injectionX, injY=injectionY)))

        # Add parameters of reservoir component model
        ares[-1].add_par('injRate', value=3.7, vary=False)
        ares[-1].add_par('logResPerm', value=-12, vary=False)
        ares[-1].add_par('reservoirRadius', value=4500.0, vary=False)
        ares[-1].add_par('brineDensity', value=1030.9, vary=False)
        ares[-1].add_par('CO2Density', value=775.0, vary=False)
        ares[-1].add_par('brineViscosity', value=7.5e-4, vary=False)
        ares[-1].add_par('CO2Viscosity', value=6.6e-5, vary=False)
        ares[-1].add_par('numberOfShaleLayers',
                         value=numberOfShaleLayers, vary=False)
        
        ares[-1].add_par_linked_to_par(
            'numberOfShaleLayers', strata[-1].deterministic_pars['numberOfShaleLayers'])
        ares[-1].add_par_linked_to_par(
            'datumPressure', strata[-1].deterministic_pars['datumPressure'])
        
        # Link the analytical reservoir's unit thicknesses to the dipping 
        # stratigraphy component's unit thicknesses.
        for ob_nm in thickness_obs:
            ares[-1].add_par_linked_to_obs(ob_nm, strata[-1].linkobs[ob_nm])
        
        ares[-1].add_obs('pressure')
        ares[-1].add_obs('CO2saturation')
        ares[-1].add_obs_to_be_linked('pressure')
        ares[-1].add_obs_to_be_linked('CO2saturation')

        # Add multisegmented wellbore component
        ms.append(sm.add_component_model_object(
            MultisegmentedWellbore(name=msName, parent=sm)))

        # Add parameters of multisegmented wellbore component model
        ms[-1].add_par('wellRadius', value=0.1, vary=False)
        ms[-1].add_par('logWellPerm', value=-11.0, vary=False)
        ms[-1].add_par('logAqu1Perm', value=-13.0, vary=False)
        ms[-1].add_par('logAqu2Perm', value=-13.0, vary=False)
        ms[-1].add_par('brineDensity', value=1030.9, vary=False)
        ms[-1].add_par('CO2Density', value=775, vary=False)
        ms[-1].add_par('brineViscosity', value=7.5e-4, vary=False)
        ms[-1].add_par('CO2Viscosity', value=6.6e-5, vary=False)

        # Add linked parameters: common to reservoir and wellbore components
        ms[-1].add_par_linked_to_par(
            'numberOfShaleLayers', strata[-1].deterministic_pars['numberOfShaleLayers'])
        ms[-1].add_par_linked_to_par(
            'datumPressure', strata[-1].deterministic_pars['datumPressure'])

        # Link the analytical reservoir's unit thicknesses to the dipping 
        # stratigraphy component's unit thicknesses.
        for ob_nm in thickness_obs:
            ms[-1].add_par_linked_to_obs(ob_nm, strata[-1].linkobs[ob_nm])

        # Add keyword arguments linked to the output provided by reservoir model
        ms[-1].add_kwarg_linked_to_obs('pressure', ares[-1].linkobs['pressure'])
        ms[-1].add_kwarg_linked_to_obs('CO2saturation', ares[-1].linkobs['CO2saturation'])

        # Add observations of multisegmented wellbore component model
        ms[-1].add_obs('brine_aquifer1')
        ms[-1].add_obs('CO2_aquifer1')
        ms[-1].add_obs('brine_aquifer2')
        ms[-1].add_obs('CO2_aquifer2')
        ms[-1].add_obs('brine_atm')
        ms[-1].add_obs('CO2_atm')

        components.append(ares[-1])
        components.append(ms[-1])

        yaml_data[ares[-1].name] = dict()
        yaml_data[ares[-1].name]['Type'] = 'AnalyticalReservoir'
        yaml_data[ms[-1].name] = dict()
        yaml_data[ms[-1].name]['Type'] = 'MultisegmentedWellbore'
        yaml_data[ms[-1].name]['Connection'] = ares[-1].name

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('pressure, Well 1\n',
          sm.collect_observations_as_time_series(ares[0], 'pressure'))
    print('------------------------------------------------------------------')
    print('pressure, Well 2\n',
          sm.collect_observations_as_time_series(ares[1], 'pressure'))
    print('------------------------------------------------------------------')
    print('CO2_aquifer1, Well 1\n',
          sm.collect_observations_as_time_series(ms[0], 'CO2_aquifer1'))
    print('------------------------------------------------------------------')
    print('brine_aquifer1, Well 2\n',
          sm.collect_observations_as_time_series(ms[1], 'brine_aquifer1'))
    print('------------------------------------------------------------------')

    if not os.path.exists(model_data['OutputDirectory']):
        os.mkdir(model_data['OutputDirectory'])

    savefig = os.path.join(model_data['OutputDirectory'], plotName)

    iam_vis.stratigraphy_plot(yaml_data, model_data, sm,
                              name=plotName, savefig=savefig)
