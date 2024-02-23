"""
This example demonstrates a system model utilizing Stratigraphy, AnalyticalReservoir, 
MultisegmentedWellboreAI, RateToMassAdapter, and GenericAquifer components.

The get_composite_depth_names() and get_depth_expr() methods of the Stratigraphy 
component are used to calculate the composite parameters related to unit depths 
(bottom, middle, and top depths for all shales and aquifers as well as the reservoir).

The pressures and CO2 saturations produced by the AnalyticalReservoir are provided 
to the MultisegmentedWellboreAI, which then produces leakage rates. The leakage 
rates produced by the MultisegmentedWellboreAI are converted into leaked masses 
with the RateToMassAdapter. Finally, these leaked masses are provided to the 
GenericAquifer, which simulates the development of dissolved salt and dissolved 
CO2 plumes in the selected aquifer. This example only focuses on leakage into 
and plume development within one aquifer. The aquifer selected can be altered 
with the selected_aquifer_num variable, which can range from 0 to 
(len(aquifer_thicknesses) - 1).

After the simulation has been run, the results are plotted.

To run this example, the files for the Multisegmented Wellbore AI component must be 
downloaded. The files are not included in the download of NRAP-Open-IAM because they are 
large (about 7.2 GB). The two files can be downloaded from these links:
    http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/42882a99-a5b0-4d68-8f4c-cad4812358c4/download
    http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/8b786e28-0e52-40d4-84b5-23c25dc2d338/download

The zipped files must then be unzipped, with the resulting folders placed into this directory:
    src/openiam/components/models/wellbore/multisegmented_ai

Example of run:
$ python iam_sys_strata_reservoir_mswellai_genericaquifer.py
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
    shale_thicknesses = [322, 216, 67, 156, 54] # shales 1 through 5
    aquifer_thicknesses = [51, 45, 142, 22]     # aquifers 1 through 4
    
    # The simulation will focus on the plumes in this aquifer
    selected_aquifer_num = 1
    
    if selected_aquifer_num < 0 or selected_aquifer_num > (len(aquifer_thicknesses) - 1):
        err_msg = ''.join([
            'The selected_aquifer_num variable was set to {}. '.format(selected_aquifer_num), 
            'This variable can only range from 0 to {}. '.format(len(aquifer_thicknesses) - 1), 
            'Due to this issue, the simulation cannot continue.'])
        raise KeyError(err_msg)
    
    # Density of the brine in the reservoir, in kg/(m^3)
    brineDensity = 1050
    
    # The x and y coordinates used for the OpenWellbore. The distance between 
    # the injection location (x = 0 m, y = 0 m) and this well location should 
    # not be larger than the reservoirRadius parameter of the AnalyticalReservoir.
    # 5 km from the injection site
    locXVal = 3535.53
    locYVal = 3535.53
    
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
    
    # Add the analytical reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(
        name='ares', parent=sm, injX=0, injY=0, locX=locXVal, locY=locYVal))
    
    ares.add_par_linked_to_par('numberOfShaleLayers', 
                               strata.deterministic_pars['numberOfShaleLayers'])
    
    ares.add_par('injRate', value=0.33, vary=False)
    ares.add_par('logResPerm', value=-13.5, vary=False)
    ares.add_par('reservoirPorosity', value=0.175, vary=False)
    ares.add_par('reservoirRadius', value=10000, vary=False)
    ares.add_par('brineResSaturation', value=0.05, vary=False)
    ares.add_par('brineDensity', value=brineDensity, vary=False)
    
    # If you set up the thicknesses to vary, you would link the AnalyticalReservoir 
    # component's thicknesses to strata.pars[parameter_name] instead.
    ares.add_par_linked_to_par('reservoirThickness', 
                               strata.deterministic_pars['reservoirThickness'])
    
    for unit_num in range(1, numberOfShaleLayers + 1):
        ares.add_par_linked_to_par('shale{}Thickness'.format(unit_num), 
                                   strata.deterministic_pars[
                                       'shale{}Thickness'.format(unit_num)])
        
        if unit_num < numberOfShaleLayers:
            ares.add_par_linked_to_par('aquifer{}Thickness'.format(unit_num), 
                                       strata.deterministic_pars[
                                           'aquifer{}Thickness'.format(unit_num)])
    
    ares.add_par_linked_to_par('datumPressure',
                               strata.default_pars['datumPressure'])
    
    # Add obervations of the analytical reservoir
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    
    # Prepare the observations to be linked to the open wellbore
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')
    
    # Add multisegmented wellbore ai component
    ms = sm.add_component_model_object(
        MultisegmentedWellboreAI(name='ms', parent=sm))
    
    ms.add_par('useDLmodel', value=1, vary=False)
    ms.add_par('wellRadius', value=0.04, vary=False)
    ms.add_par('logWell1Perm', value=-12.5, vary=False) # perm. for shale 1
    ms.add_par('logAqu1Perm', value=-13.5, vary=False)  # perm. for aquifer 1
    ms.add_par('logWell2Perm', value=-15.5, vary=False) # perm. for shale 2
    ms.add_par('logAqu2Perm', value=-13.5, vary=False)  # perm. for aquifer 2
    ms.add_par('logWell3Perm', value=-15.5, vary=False) # perm. for shale 3
    ms.add_par('logAqu3Perm', value=-13.5, vary=False)  # perm. for aquifer 3
    ms.add_par('logWell4Perm', value=-15.5, vary=False) # perm. for shale 4
    ms.add_par('logAqu4Perm', value=-13.5, vary=False)  # perm. for aquifer 4
    ms.add_par('logWell5Perm', value=-15.5, vary=False) # perm. for shale 5
    
    # Add linked parameters: common to stratigraphy and wellbore components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             strata.deterministic_pars['numberOfShaleLayers'])
    
    for unit_num in range(1, numberOfShaleLayers + 1):
        ms.add_par_linked_to_par('shale{}Thickness'.format(unit_num), 
                                 strata.deterministic_pars[
                                     'shale{}Thickness'.format(unit_num)])
        
        if unit_num < numberOfShaleLayers:
            ms.add_par_linked_to_par('aquifer{}Thickness'.format(unit_num), 
                                     strata.deterministic_pars[
                                         'aquifer{}Thickness'.format(unit_num)])
    
    ms.add_par_linked_to_par('reservoirThickness', 
                             strata.deterministic_pars['reservoirThickness'])
    
    ms.add_par_linked_to_par('datumPressure', strata.default_pars['datumPressure'])
    
    # Add linked parameters: common to stratigraphy and wellbore components
    ms.add_par_linked_to_par('brineDensity', ares.deterministic_pars['brineDensity'])
    
    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])
    
    # Add observations of multisegmented wellbore component AI model
    for unit_num in range(1, numberOfShaleLayers):
        ms.add_obs_to_be_linked('brine_aquifer{}'.format(unit_num))
        ms.add_obs_to_be_linked('CO2_aquifer{}'.format(unit_num))
        
        ms.add_obs('brine_aquifer{}'.format(unit_num))
        ms.add_obs('CO2_aquifer{}'.format(unit_num))
    
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs_to_be_linked('CO2_atm')
    
    ms.add_obs('brine_atm')
    ms.add_obs('CO2_atm')
    
    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(
        RateToMassAdapter(name='adapt', parent=sm))
    
    adapt.add_kwarg_linked_to_collection('brine_aquifer1',
        [ms.linkobs['brine_aquifer1'], ms.linkobs['brine_aquifer2']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer2',
        [ms.linkobs['brine_aquifer2'], ms.linkobs['brine_aquifer3']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer3',
        [ms.linkobs['brine_aquifer3'], ms.linkobs['brine_aquifer4']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer4',
        [ms.linkobs['brine_aquifer4'], ms.linkobs['brine_atm']])
    
    adapt.add_kwarg_linked_to_collection('CO2_aquifer1',
        [ms.linkobs['CO2_aquifer1'], ms.linkobs['CO2_aquifer2']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer2',
        [ms.linkobs['CO2_aquifer2'], ms.linkobs['CO2_aquifer3']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer3',
        [ms.linkobs['CO2_aquifer3'], ms.linkobs['CO2_aquifer4']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer4',
        [ms.linkobs['CO2_aquifer4'], ms.linkobs['CO2_atm']])
    
    # Add observations of RateToMassAdapter model
    for unit_num in range(1, numberOfShaleLayers):
        adapt.add_obs_to_be_linked('mass_brine_aquifer{}'.format(unit_num))
        adapt.add_obs_to_be_linked('mass_CO2_aquifer{}'.format(unit_num))
        
        adapt.add_obs('mass_brine_aquifer{}'.format(unit_num))
        adapt.add_obs('mass_CO2_aquifer{}'.format(unit_num))
    
    # Add the generic aquifer component
    genaq = sm.add_component_model_object(GenericAquifer(name='genaq', parent=sm))
    
    # Add the parameters of the generic aquifer
    genaq.add_par('por', value=0.125, vary=False)
    genaq.add_par('log_permh', value=-13.5, vary=False)
    genaq.add_par('log_aniso', value=0.3, vary=False)
    genaq.add_par('reservoir_salinity', value=0.035, vary=False)
    genaq.add_par('aquifer_salinity', value=0.005, vary=False)
    genaq.add_par('dissolved_salt_threshold', value=0.01, vary=False)
    genaq.add_par('dissolved_co2_threshold', value=0.01, vary=False)
    
    genaq.add_par_linked_to_par('aqu_thick', strata.deterministic_pars[
        'aquifer{}Thickness'.format(selected_aquifer_num)])
    genaq.add_par_linked_to_par('top_depth', strata.composite_pars[
        'aquifer{}TopDepth'.format(selected_aquifer_num)])
    
    # Add keyword arguments of the generic aquifer component model
    genaq.add_kwarg_linked_to_obs(
        'brine_mass', adapt.linkobs['mass_brine_aquifer{}'.format(selected_aquifer_num)])
    genaq.add_kwarg_linked_to_obs(
        'co2_mass', adapt.linkobs['mass_CO2_aquifer{}'.format(selected_aquifer_num)])
    
    # Add observations of the generic aquifer
    genaq.add_obs('Dissolved_salt_volume')
    genaq.add_obs('Dissolved_CO2_volume')
    genaq.add_obs('Dissolved_salt_dr')
    genaq.add_obs('Dissolved_CO2_dr')
    genaq.add_obs('Dissolved_salt_dz')
    genaq.add_obs('Dissolved_CO2_dz')
    
    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    
    # Collect observations
    pressure = sm.collect_observations_as_time_series(ares, 'pressure')
    CO2saturation = sm.collect_observations_as_time_series(ares, 'CO2saturation')
    
    brine_aquifer1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')
    CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
    
    brine_aquifer2 = sm.collect_observations_as_time_series(ms, 'brine_aquifer2')
    CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer2')
    
    brine_aquifer3 = sm.collect_observations_as_time_series(ms, 'brine_aquifer3')
    CO2_aquifer3 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer3')
    
    brine_aquifer4 = sm.collect_observations_as_time_series(ms, 'brine_aquifer4')
    CO2_aquifer4 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer4')
    
    mass_brine_aquifer1 = sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer1')
    mass_CO2_aquifer1 = sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer1')
    
    mass_brine_aquifer2 = sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer2')
    mass_CO2_aquifer2 = sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer2')
    
    mass_brine_aquifer3 = sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer3')
    mass_CO2_aquifer3 = sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer3')
    
    mass_brine_aquifer4 = sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer4')
    mass_CO2_aquifer4 = sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer4')
    
    Dissolved_salt_volume = sm.collect_observations_as_time_series(
        genaq, 'Dissolved_salt_volume')
    Dissolved_CO2_volume = sm.collect_observations_as_time_series(
        genaq, 'Dissolved_CO2_volume')
    
    Dissolved_salt_dr = sm.collect_observations_as_time_series(
        genaq, 'Dissolved_salt_dr')
    Dissolved_CO2_dr = sm.collect_observations_as_time_series(
        genaq, 'Dissolved_CO2_dr')
    
    Dissolved_salt_dz = sm.collect_observations_as_time_series(
        genaq, 'Dissolved_salt_dz')
    Dissolved_CO2_dz = sm.collect_observations_as_time_series(
        genaq, 'Dissolved_CO2_dz')
    
    # Make the figures
    genfontsize = 10
    lgndfontsize = 10
    axisfontsize = 12
    titlefontsize = 14
    
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': genfontsize}
    plt.rc('font', **font)
    
    figsize = (14, 8)
    grid_alpha_val = 0.25
    dpi_ref = 300
    linewidth_ref = 3
    time_array_years = time_array / 365.25
    
    # Reservoir Conditions
    plt.figure(1, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Reservoir Conditions Over Time', fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(1,2,1)
    ax.plot(time_array_years, pressure, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Pressure (Pa)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    
    ax = plt.subplot(1,2,2)
    ax.plot(time_array_years, CO2saturation, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('CO$_2$ Saturation', fontsize=axisfontsize, fontweight='bold')
    
    # Brine Leakage Rates
    plt.figure(2, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Brine Leakage Rates', fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,1)
    ax.plot(time_array_years, brine_aquifer1, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_ylabel('Brine Leakage Rate (kg/s)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 1', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,2)
    ax.plot(time_array_years, brine_aquifer2, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 2', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,3)
    ax.plot(time_array_years, brine_aquifer3, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Brine Leakage Rate (kg/s)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 3', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,4)
    ax.plot(time_array_years, brine_aquifer4, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Brine Leakage Rate (kg/s)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 4', fontsize=axisfontsize, fontweight='bold')
    
    # CO2 Leakage Rates
    plt.figure(3, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('CO$_2$ Leakage Rates', fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,1)
    ax.plot(time_array_years, CO2_aquifer1, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_ylabel('CO$_2$ Leakage Rate (kg/s)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 1', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,2)
    ax.plot(time_array_years, CO2_aquifer2, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 2', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,3)
    ax.plot(time_array_years, CO2_aquifer3, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('CO$_2$ Leakage Rate (kg/s)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 3', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,4)
    ax.plot(time_array_years, CO2_aquifer4, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 4', fontsize=axisfontsize, fontweight='bold')
    
    # Leaked Brine Masses
    plt.figure(4, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Cumulative Brine Leakage', fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,1)
    ax.plot(time_array_years, mass_brine_aquifer1, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 1', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,2)
    ax.plot(time_array_years, mass_brine_aquifer2, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 2', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,3)
    ax.plot(time_array_years, mass_brine_aquifer3, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Leaked Brine Mass (kg)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 3', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,4)
    ax.plot(time_array_years, mass_brine_aquifer4, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 4', fontsize=axisfontsize, fontweight='bold')
    
    # Leaked CO2 Masses
    plt.figure(5, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Cumulative CO$_2$ Leakage', fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,1)
    ax.plot(time_array_years, mass_CO2_aquifer1, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_ylabel('Leaked CO$_2$ Mass (kg)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 1', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,2)
    ax.plot(time_array_years, mass_CO2_aquifer2, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 2', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,3)
    ax.plot(time_array_years, mass_CO2_aquifer3, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Leaked CO$_2$ Mass (kg)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 3', fontsize=axisfontsize, fontweight='bold')
    
    ax = plt.subplot(2,2,4)
    ax.plot(time_array_years, mass_CO2_aquifer4, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText='True')
    ax.set_title('Aquifer 4', fontsize=axisfontsize, fontweight='bold')
    
    # Plume Volumes
    plt.figure(6, figsize=figsize, dpi=dpi_ref)
    
    plt.suptitle('Impacts on Aquifer {}'.format(selected_aquifer_num), 
                 fontsize=titlefontsize, fontweight='bold')
    
    ax = plt.subplot(1,2,1)
    
    ax.plot(time_array_years, Dissolved_salt_volume, linewidth=linewidth_ref)
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Dissolved Salt Volume (m$^3$)', fontsize=axisfontsize, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    
    ax = plt.subplot(1,2,2)
    
    ax.plot(time_array_years, Dissolved_CO2_volume, linewidth=linewidth_ref)
    
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
    
    ax = plt.subplot(1,2,1)
    
    ax.plot(time_array_years, Dissolved_salt_dr, linewidth=linewidth_ref, 
            color='C0', label='Radius')
    ax.plot(time_array_years, Dissolved_salt_dz, linewidth=linewidth_ref, 
            color='C1', linestyle=':', label='Height')
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Dissolved Salt Plume Radius or Height (m)', fontsize=axisfontsize, 
                  fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    plt.legend(fontsize=lgndfontsize)
    
    ax = plt.subplot(1,2,2)
    
    ax.plot(time_array_years, Dissolved_CO2_dr, linewidth=linewidth_ref, 
            color='C0', label='Radius')
    ax.plot(time_array_years, Dissolved_CO2_dz, linewidth=linewidth_ref, 
            color='C1', linestyle=':', label='Height')
    
    ax.grid(alpha=grid_alpha_val)
    ax.set_xlabel('Time (years)', fontsize=axisfontsize, fontweight='bold')
    ax.set_ylabel('Dissolved CO$_2$ Plume Radius or Height (m)', fontsize=axisfontsize, 
                  fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText='True')
    plt.legend(fontsize=lgndfontsize)
