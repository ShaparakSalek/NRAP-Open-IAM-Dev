"""
This example is similar to the example iam_sys_strata_reservoir_mswellai_genericaquifer.py, 
but this example evaluates leakged through multiple MultisegmentedWellboreAI components. The 
x and y coordinates of these components are defined in the lists locXVal and 
locYVal. The AnalyticalReservoir, OpenWellbore, RateToMassAdapter, and GenericAquifer 
components made for each location are stored in lists.

After the simulation has been run, the results are stored in lists and plotted.

Because this model uses 5 wellbore locations and using the machine learning models is computationally 
demanding, the useDLmodel parameter is set to 0. The Multisegmented Wellbore AI results are therefore 
derived from analytical modeling.

To run this example, one of the files for the Multisegmented Wellbore AI component must be 
downloaded. The files are not included in the download of NRAP-Open-IAM because they are 
large. The two files can be downloaded from these links:
    http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/42882a99-a5b0-4d68-8f4c-cad4812358c4/download
    http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/8b786e28-0e52-40d4-84b5-23c25dc2d338/download

Because useDLmodel is set to 0 here, only the file from the first link is required (about 2 GB). Although  
the results are based on analytical modeling in this case, the file is still required to calculate 
fluid properties. The second file (about 5.2 GB) is required if useDLmodel is set to 1.

The zipped file must then be unzipped, with the resulting folder placed into this directory:
    src/openiam/components/models/wellbore/multisegmented_ai

Example of run:
$ python iam_sys_strata_reservoir_mswellai_genericaquifer_5locs.py
"""

import logging
import numpy as np
import matplotlib.pyplot as plt

from openiam.components.iam_base_classes import SystemModel
from openiam.components.stratigraphy_component import Stratigraphy
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.multisegmented_wellbore_ai_component import MultisegmentedWellboreAI
from openiam.components.rate_to_mass_adapter import RateToMassAdapter
from openiam.components.generic_aquifer_component import GenericAquifer


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    # Define keyword arguments of the system model. The times used do not have 
    # an impact on the composite depth parameters printed further below.
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    
    numberOfShaleLayers = 5
    reservoir_thickness = 31
    shale_thicknesses = [322, 216, 67, 156, 54]  # shales 1 through 5
    aquifer_thicknesses = [51, 45, 142, 22]      # aquifers 1 through 4
    
    # The simulation will focus on the brine and CO2 leakage into this aquifer
    selected_aquifer_num = 1
    
    if selected_aquifer_num < 0 or selected_aquifer_num > (len(aquifer_thicknesses) - 1):
        err_msg = ''.join([
            'The selected_aquifer_num variable was set to {}. '.format(selected_aquifer_num), 
            'This variable can only range from 0 to {}. '.format(len(aquifer_thicknesses) - 1), 
            'Due to this issue, the simulation cannot continue.'])
        raise KeyError(err_msg)
    
    # Density of the brine in the reservoir, in kg/(m^3)
    brineDensity = 1050
    
    # The x and y coordinates used for the OpenWellbores. The distance between 
    # the injection location (x = 0 m, y = 0 m) and any of these well locations 
    # should not be larger than the reservoirRadius parameter of the AnalyticalReservoir.
    # Distances of 5 km to 9 km from x = 0 m, y = 0 m.
    locXVal = [3535.53, 4242.64, 4949.75, 5656.85, 6363.96]
    locYVal = [3535.53, 4242.64, 4949.75, 5656.85, 6363.96]
    
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
    
    # Lists containing the components
    ares = []
    mws = []
    adapts = []
    genaqs = []
    # Make the components for each location
    for locRef in range(len(locXVal)):
        # Add the analytical reservoir component
        ares.append(sm.add_component_model_object(AnalyticalReservoir(
            name='ares{}'.format(locRef), parent=sm, injX=0, injY=0, 
            locX=locXVal[locRef], locY=locYVal[locRef])))
        
        ares[-1].add_par_linked_to_par(
            'numberOfShaleLayers', strata.deterministic_pars['numberOfShaleLayers'])
        
        ares[-1].add_par('injRate', value=0.5, vary=False)
        ares[-1].add_par('logResPerm', value=-13.5, vary=False)
        ares[-1].add_par('reservoirPorosity', value=0.175, vary=False)
        ares[-1].add_par('reservoirRadius', value=10000, vary=False)
        ares[-1].add_par('brineResSaturation', value=0.05, vary=False)
        ares[-1].add_par('brineDensity', value=brineDensity, vary=False)
        
        # If you set up the thicknesses to vary, you would link the AnalyticalReservoir 
        # component's thicknesses to strata.pars[parameter_name] instead.
        ares[-1].add_par_linked_to_par(
            'reservoirThickness', strata.deterministic_pars['reservoirThickness'])
        
        for unit_num in range(1, numberOfShaleLayers + 1):
            ares[-1].add_par_linked_to_par(
                'shale{}Thickness'.format(unit_num), strata.deterministic_pars[
                    'shale{}Thickness'.format(unit_num)])
            
            if unit_num < numberOfShaleLayers:
                ares[-1].add_par_linked_to_par(
                    'aquifer{}Thickness'.format(unit_num), strata.deterministic_pars[
                        'aquifer{}Thickness'.format(unit_num)])
        
        ares[-1].add_par_linked_to_par('datumPressure', 
                                       strata.default_pars['datumPressure'])
        
        # Add obervations of the analytical reservoir
        ares[-1].add_obs('pressure')
        ares[-1].add_obs('CO2saturation')
        
        # Prepare the observations to be linked to the open wellbore
        ares[-1].add_obs_to_be_linked('pressure')
        ares[-1].add_obs_to_be_linked('CO2saturation')
        
        # Add the multisegmented wellbore ai component
        mws.append(sm.add_component_model_object(
            MultisegmentedWellboreAI(name='mswai{}'.format(locRef), parent=sm)))
        
        mws[-1].add_par('useDLmodel', value=0, vary=False)
        mws[-1].add_par('wellRadius', value=0.04, vary=False)
        mws[-1].add_par('logWell1Perm', value=-12.5, vary=False) # shale 1
        mws[-1].add_par('logAqu1Perm', value=-13.0, vary=False)  # aquifer 1
        mws[-1].add_par('logWell2Perm', value=-15.5, vary=False) # shale 2
        mws[-1].add_par('logAqu2Perm', value=-13.5, vary=False)  # aquifer 2
        mws[-1].add_par('logWell3Perm', value=-15.5, vary=False) # shale 3
        mws[-1].add_par('logAqu3Perm', value=-13.5, vary=False)  # aquifer 3
        mws[-1].add_par('logWell4Perm', value=-15.5, vary=False) # shale 4
        mws[-1].add_par('logAqu4Perm', value=-13.5, vary=False)  # aquifer 4
        mws[-1].add_par('logWell5Perm', value=-15.5, vary=False) # shale 5
        
        # Add linked parameters: common to stratigraphy and wellbore components
        mws[-1].add_par_linked_to_par(
            'numberOfShaleLayers', strata.deterministic_pars['numberOfShaleLayers'])
        
        for unit_num in range(1, numberOfShaleLayers + 1):
            mws[-1].add_par_linked_to_par(
                'shale{}Thickness'.format(unit_num), strata.deterministic_pars[
                    'shale{}Thickness'.format(unit_num)])
            
            if unit_num < numberOfShaleLayers:
                mws[-1].add_par_linked_to_par(
                    'aquifer{}Thickness'.format(unit_num), strata.deterministic_pars[
                        'aquifer{}Thickness'.format(unit_num)])
        
        mws[-1].add_par_linked_to_par(
            'reservoirThickness', strata.deterministic_pars['reservoirThickness'])
        
        mws[-1].add_par_linked_to_par(
            'datumPressure', strata.default_pars['datumPressure'])

        # Add linked parameters: common to stratigraphy and wellbore components
        mws[-1].add_par_linked_to_par(
            'brineDensity', ares[-1].deterministic_pars['brineDensity'])
        
        # Add keyword arguments linked to the output provided by reservoir model
        mws[-1].add_kwarg_linked_to_obs(
            'pressure', ares[-1].linkobs['pressure'])
        mws[-1].add_kwarg_linked_to_obs(
            'CO2saturation', ares[-1].linkobs['CO2saturation'])
        
        # Add observations of multisegmented wellbore component AI model
        for unit_num in range(1, numberOfShaleLayers):
            mws[-1].add_obs_to_be_linked('brine_aquifer{}'.format(unit_num))
            mws[-1].add_obs_to_be_linked('CO2_aquifer{}'.format(unit_num))
            
            mws[-1].add_obs('brine_aquifer{}'.format(unit_num))
            mws[-1].add_obs('CO2_aquifer{}'.format(unit_num))

        mws[-1].add_obs_to_be_linked('brine_atm')
        mws[-1].add_obs_to_be_linked('CO2_atm')
        
        mws[-1].add_obs('brine_atm')
        mws[-1].add_obs('CO2_atm')

        # Add adapter that transforms leakage rates to accumulated mass
        adapts.append(sm.add_component_model_object(
            RateToMassAdapter(name='adapt{}'.format(locRef), parent=sm)))
        
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer1',
            [mws[-1].linkobs['CO2_aquifer1'], mws[-1].linkobs['CO2_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer2',
            [mws[-1].linkobs['CO2_aquifer2'], mws[-1].linkobs['CO2_aquifer3']])
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer3',
            [mws[-1].linkobs['CO2_aquifer3'], mws[-1].linkobs['CO2_aquifer4']])
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer4',
            [mws[-1].linkobs['CO2_aquifer4'], mws[-1].linkobs['CO2_atm']])
        
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer1',
            [mws[-1].linkobs['brine_aquifer1'], mws[-1].linkobs['brine_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer2',
            [mws[-1].linkobs['brine_aquifer2'], mws[-1].linkobs['brine_aquifer3']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer3',
            [mws[-1].linkobs['brine_aquifer3'], mws[-1].linkobs['brine_aquifer4']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer4',
            [mws[-1].linkobs['brine_aquifer4'], mws[-1].linkobs['brine_atm']])
        
        # Add observations of RateToMassAdapter model
        for unit_num in range(1, numberOfShaleLayers):
            adapts[-1].add_obs_to_be_linked('mass_brine_aquifer{}'.format(unit_num))
            adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer{}'.format(unit_num))
            
            adapts[-1].add_obs('mass_brine_aquifer{}'.format(unit_num))
            adapts[-1].add_obs('mass_CO2_aquifer{}'.format(unit_num))
        
        # Add the generic aquifer component
        genaqs.append(sm.add_component_model_object(GenericAquifer(
            name='genaq{}'.format(locRef), parent=sm)))
        
        # Add the parameters of the generic aquifer
        genaqs[-1].add_par('por', value=0.125, vary=False)
        genaqs[-1].add_par('log_permh', value=-13.0, vary=False)
        genaqs[-1].add_par('log_aniso', value=0.3, vary=False)
        genaqs[-1].add_par('reservoir_salinity', value=0.035, vary=False)
        genaqs[-1].add_par('aquifer_salinity', value=0.005, vary=False)
        genaqs[-1].add_par('dissolved_salt_threshold', value=0.01, vary=False)
        genaqs[-1].add_par('dissolved_co2_threshold', value=0.01, vary=False)
        
        genaqs[-1].add_par_linked_to_par('aqu_thick', strata.deterministic_pars[
            'aquifer{}Thickness'.format(selected_aquifer_num)])
        genaqs[-1].add_par_linked_to_par('top_depth', strata.composite_pars[
            'aquifer{}TopDepth'.format(selected_aquifer_num)])
        
        # Add keyword arguments of the generic aquifer component model
        genaqs[-1].add_kwarg_linked_to_obs(
            'brine_mass', adapts[-1].linkobs['mass_brine_aquifer{}'.format(selected_aquifer_num)])
        genaqs[-1].add_kwarg_linked_to_obs(
            'co2_mass', adapts[-1].linkobs['mass_CO2_aquifer{}'.format(selected_aquifer_num)])
        
        # Add observations of the generic aquifer
        genaqs[-1].add_obs('Dissolved_salt_volume')
        genaqs[-1].add_obs('Dissolved_CO2_volume')
        genaqs[-1].add_obs('Dissolved_salt_dr')
        genaqs[-1].add_obs('Dissolved_CO2_dr')
        genaqs[-1].add_obs('Dissolved_salt_dz')
        genaqs[-1].add_obs('Dissolved_CO2_dz')
    
    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    
    # Collect the observations in lists
    pressure = []
    CO2saturation = []
    brine_aquifers = []
    CO2_aquifers = []
    mass_brine_aquifers = []
    mass_CO2_aquifers = []
    Dissolved_salt_volume = []
    Dissolved_CO2_volume = []
    Dissolved_salt_dr = []
    Dissolved_CO2_dr = []
    Dissolved_salt_dz = []
    Dissolved_CO2_dz = []
    
    for locRef in range(len(locXVal)):
        pressure.append(sm.collect_observations_as_time_series(
            ares[locRef], 'pressure'))
        CO2saturation.append(sm.collect_observations_as_time_series(
            ares[locRef], 'CO2saturation'))
        
        brine_aquifers.append([
            sm.collect_observations_as_time_series(
                mws[locRef], 'brine_aquifer{}'.format(aquifer_num)) 
            for aquifer_num in range(1, numberOfShaleLayers)])
        
        CO2_aquifers.append([
            sm.collect_observations_as_time_series(
                mws[locRef], 'CO2_aquifer{}'.format(aquifer_num)) 
            for aquifer_num in range(1, numberOfShaleLayers)])
        
        mass_brine_aquifers.append([
            sm.collect_observations_as_time_series(
                adapts[locRef], 'mass_brine_aquifer{}'.format(aquifer_num)) 
            for aquifer_num in range(1, numberOfShaleLayers)])
        
        mass_CO2_aquifers.append([
            sm.collect_observations_as_time_series(
                adapts[locRef], 'mass_CO2_aquifer{}'.format(aquifer_num)) 
            for aquifer_num in range(1, numberOfShaleLayers)])
        
        Dissolved_salt_volume.append(sm.collect_observations_as_time_series(
            genaqs[locRef], 'Dissolved_salt_volume'))
        Dissolved_CO2_volume.append(sm.collect_observations_as_time_series(
            genaqs[locRef], 'Dissolved_CO2_volume'))
        
        Dissolved_salt_dr.append(sm.collect_observations_as_time_series(
            genaqs[locRef], 'Dissolved_salt_dr'))
        Dissolved_CO2_dr.append(sm.collect_observations_as_time_series(
            genaqs[locRef], 'Dissolved_CO2_dr'))
        
        Dissolved_salt_dz.append(sm.collect_observations_as_time_series(
            genaqs[locRef], 'Dissolved_salt_dz'))
        Dissolved_CO2_dz.append(sm.collect_observations_as_time_series(
            genaqs[locRef], 'Dissolved_CO2_dz'))
    
    # Make the figures
    genfontsize = 10
    lgndfontsize = 10
    axisfontsize = 12
    titlefontsize = 14
    
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': genfontsize}
    plt.rc('font', **font)
    
    figsize = (14,8)
    grid_alpha_val = 0.25
    dpi_ref = 300
    linewidth_ref = 3
    time_array_years = time_array / 365.25
    
    subplots_with_xlabel = [3, 4]
    subplots_with_ylabel = [1, 3]
    
    # Reservoir Conditions
    plt.figure(1, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Reservoir Conditions Over Time', fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(1,2,1)
    
    for locRef in range(len(locXVal)):
        ax.plot(time_array_years, pressure[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
    
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Pressure (Pa)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    
    ax = plt.subplot(1,2,2)
    for locRef in range(len(locXVal)):
        ax.plot(time_array_years, CO2saturation[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
    
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('CO$_2$ Saturation', fontsize=axisfontsize, fontweight='bold')
    
    # Brine Leakage Rates
    plt.figure(2, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Brine Leakage Rates', fontsize=titlefontsize, fontweight='bold')
    
    for aquifer_num in range(1, numberOfShaleLayers):
        ax = plt.subplot(2,2,aquifer_num)
        
        for locRef in range(len(locXVal)):
            ax.plot(time_array_years, brine_aquifers[locRef][aquifer_num - 1], linewidth=linewidth_ref, 
                    color='C{}'.format(locRef), label='Location {}'.format(locRef))
        
        ax.legend(fontsize=lgndfontsize)
        ax.grid(alpha=grid_alpha_val)
        if aquifer_num in subplots_with_xlabel:
            ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
        if aquifer_num in subplots_with_ylabel:
            ax.set_ylabel('Brine Leakage Rate (kg/s)', fontsize=axisfontsize, fontweight='bold')
        ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
        ax.set_title('Aquifer {}'.format(aquifer_num), fontsize=axisfontsize, fontweight='bold')
    
    # CO2 Leakage Rates
    plt.figure(3, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('CO$_2$ Leakage Rates', fontsize=titlefontsize, fontweight='bold')
    
    for aquifer_num in range(1, numberOfShaleLayers):
        ax = plt.subplot(2,2,aquifer_num)
        
        for locRef in range(len(locXVal)):
            ax.plot(time_array_years, CO2_aquifers[locRef][aquifer_num - 1], linewidth=linewidth_ref, 
                    color='C{}'.format(locRef), label='Location {}'.format(locRef))
        
        ax.legend(fontsize=lgndfontsize)
        ax.grid(alpha=grid_alpha_val)
        if aquifer_num in subplots_with_xlabel:
            ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
        if aquifer_num in subplots_with_ylabel:
            ax.set_ylabel('CO$_2$ Leakage Rate (kg/s)', fontsize=axisfontsize, fontweight='bold')
        ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
        ax.set_title('Aquifer {}'.format(aquifer_num), fontsize=axisfontsize, fontweight='bold')
    
    # Leaked Brine Masses
    plt.figure(4, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Cumulative Brine Leakage', fontsize=titlefontsize, fontweight='bold')
    
    for aquifer_num in range(1, numberOfShaleLayers):
        ax = plt.subplot(2,2,aquifer_num)
        
        for locRef in range(len(locXVal)):
            ax.plot(time_array_years, mass_brine_aquifers[locRef][aquifer_num - 1], linewidth=linewidth_ref, 
                    color='C{}'.format(locRef), label='Location {}'.format(locRef))
        
        ax.legend(fontsize=lgndfontsize)
        ax.grid(alpha=grid_alpha_val)
        if aquifer_num in subplots_with_xlabel:
            ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
        if aquifer_num in subplots_with_ylabel:
            ax.set_ylabel('Leaked Brine Mass (kg)', fontsize=axisfontsize, fontweight='bold')
        ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
        ax.set_title('Aquifer {}'.format(aquifer_num), fontsize=axisfontsize, fontweight='bold')
    
    # Leaked CO2 Masses
    plt.figure(5, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Cumulative CO$_2$ Leakage', fontsize=titlefontsize, fontweight='bold')
    
    for aquifer_num in range(1, numberOfShaleLayers):
        ax = plt.subplot(2,2,aquifer_num)
        
        for locRef in range(len(locXVal)):
            ax.plot(time_array_years, mass_CO2_aquifers[locRef][aquifer_num - 1], linewidth=linewidth_ref, 
                    color='C{}'.format(locRef), label='Location {}'.format(locRef))
        
        ax.legend(fontsize=lgndfontsize)
        ax.grid(alpha=grid_alpha_val)
        if aquifer_num in subplots_with_xlabel:
            ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
        if aquifer_num in subplots_with_ylabel:
            ax.set_ylabel('Leaked CO$_2$ Mass (kg)', fontsize=axisfontsize, fontweight='bold')
        ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
        ax.set_title('Aquifer {}'.format(aquifer_num), fontsize=axisfontsize, fontweight='bold')
    
    # Plume Volumes
    plt.figure(6, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Impacts on Aquifer {}'.format(selected_aquifer_num), 
                  fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(1,2,1)
    for locRef in range(len(locXVal)):
        ax.plot(time_array_years, Dissolved_salt_volume[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
    
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Dissolved Salt Volume (m$^3$)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    
    ax = plt.subplot(1,2,2)
    for locRef in range(len(locXVal)):
        ax.plot(time_array_years, Dissolved_CO2_volume[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
    
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Dissolved CO$_2$ Volume (m$^3$)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    
    # Plume Dimensions
    plt.figure(7, figsize=figsize, dpi=dpi_ref)
    
    selected_aquifer_thickness = strata.deterministic_pars[
        'aquifer{}Thickness'.format(selected_aquifer_num)].value
    
    plt.suptitle('Plume Dimensions Over Time,\nAquifer {} with a Thickness of {} m'.format(
        selected_aquifer_num, selected_aquifer_thickness), fontsize=titlefontsize, 
        fontweight='bold')
    
    for locRef in range(len(locXVal)):
        ax = plt.subplot(2,2,1)
        ax.plot(time_array_years, Dissolved_salt_dr[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
        
        ax = plt.subplot(2,2,2)
        ax.plot(time_array_years, Dissolved_salt_dz[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
    
    ax = plt.subplot(2,2,1)
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_ylabel('Dissolved Salt Plume Radius (m)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
    
    ax = plt.subplot(2,2,2)
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_ylabel('Dissolved Salt Plume Height (m)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
    
    for locRef in range(len(locXVal)):
        ax = plt.subplot(2,2,3)
        ax.plot(time_array_years, Dissolved_CO2_dr[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
        
        ax = plt.subplot(2,2,4)
        ax.plot(time_array_years, Dissolved_CO2_dz[locRef], linewidth=linewidth_ref, 
                color='C{}'.format(locRef), label='Location {}'.format(locRef))
    
    ax = plt.subplot(2,2,3)
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Dissolved CO$_2$ Plume Radius (m)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
    
    ax = plt.subplot(2,2,4)
    ax.legend(fontsize=lgndfontsize)
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Dissolved CO$_2$ Plume Height (m)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',scilimits=(0, 0), useMathText='True')
