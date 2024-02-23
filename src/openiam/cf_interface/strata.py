"""
Created: September 9th, 2022
Last modified: February 21st, 2024

Authors: Nate Mitchell, Veronika Vasylkivska
Leidos supporting NETL
"""

import logging
import numpy as np
from matplotlib.colors import is_color_like

from openiam.components.stratigraphy_component import Stratigraphy
from openiam.components.dipping_stratigraphy_component import DippingStratigraphy
from openiam.components.lookup_table_stratigraphy_component import LookupTableStratigraphy
import openiam.cf_interface.commons as cfi_commons

# Update this dictionary if the components are changed or new stratigraphy components are added 
DEFAULT_NUM_SHALE_LAYERS_DICT = {
    'Stratigraphy': 3, 'DippingStratigraphy': 3, 'LookupTableStratigraphy': 3}

DEFAULT_NUM_SHALE_LAYERS = 3

# Functions in this file are used to generate the colors, alpha values, and
# labels for different units.
# Color used to plot the unit
UNIT_COLOR_DICT = {
    'ReservoirColor': [0.5, 0.5, 0.5],
    'ShaleColor': [1, 0, 0],
    'AquiferColor': [0, 0, 1],
    }

# Default well color and alpha values
WELL_COLOR = [0, 0, 1]
WELL_ALPHA = 1
WELL_ALPHA_FILL = 0.75

# This is used to scale how much higher the alpha of labels / lines is than
# the alpha for filled areas:
# alpha = alphaFill + ((1 - alphaFill) * ALPHA_SCALE_FACTOR)
ALPHA_SCALE_FACTOR = 0.8

# Alpha used when plotting unit labels and the lines on the edges of units
UNIT_ALPHA_DICT = {
    'ReservoirAlpha': 0.9,
    'ShaleAlpha': 0.85,
    'AquiferAlpha': 0.85,
    }

# Alpha used when plotting a filled in area
UNIT_ALPHA_FILL_DICT = {
    'ReservoirAlphaFill': 0.5,
    'ShaleAlphaFill': 0.25,
    'AquiferAlphaFill': 0.25,
    }

# Labels for each unit. This is used in stratigraphy_plot.py (thickness labels
# are too long for that plot type).
UNIT_LABEL_DICT = {
    'ReservoirLabel': 'Reservoir',
    'ShaleLabel': 'Shale {:.0f}',
    'AquiferLabel': 'Aquifer {:.0f}',
    }

# Labels that include the unit thickness. This is used in stratigraphic_column.py.
UNIT_LABEL_THICKNESS_DICT = {
    'ReservoirLabel': 'Reservoir Thickness: {:.2f} m',
    'ShaleLabel': 'Shale {:.0f} Thickness: {:.2f} m',
    'AquiferLabel': 'Aquifer {:.0f} Thickness: {:.2f} m',
    }

DIP_DIRECTION_DEGREE_OPTIONS = [0, 90, 180, 270, 360]
DIP_DIRECTION_OPTIONS = ['N', 'S', 'E', 'W', 'NE', 'NW', 'SE', 'SW']

CHECK_CONDITIONS_MSG = ''.join([
    'Check your input. Alternatively, calculating unit thicknesses with these ',
    'strike and dip values may not be appropriate for these locations and/or ',
    'conditions.'])


def get_comp_types_strata_pars():
    """
    Returns a list of the stratigraphy component types that offer thickness and
    depth values as parameters (including composite parameters).
    """
    comp_types = ['Stratigraphy']

    return comp_types


def get_comp_types_strata_obs():
    """
    Returns a list of the stratigraphy component types that offer thickness and
    depth values as observations.
    """
    comp_types = ['LookupTableStratigraphy', 'DippingStratigraphy']

    return comp_types


def initialize_strata(yaml_data, sm):
    """
    This function handles the processing of stratigraphy input from the .yaml
    file.

    :param yaml_data: dictionary of input values from the .yaml file
    :type yaml_data: dict

    :param sm: system model object
    :type sm: openiam.iam_base_classes.SystemModel

    :returns: strata, sm, strata_type
    """
    # Indicates the component type
    strata_type = get_strata_type_from_yaml(yaml_data)

    types_strata_pars = get_comp_types_strata_pars()

    strata = []

    if strata_type in types_strata_pars:
        # Add Stratigraphy component
        strata = sm.add_component_model_object(Stratigraphy(
            name='strata', parent=sm))

        comp_data = yaml_data['Stratigraphy']

        if 'Parameters' in comp_data:
            strat_params = comp_data['Parameters']
        else:
            strat_params = comp_data

        # Add the parameters
        for parameter in strat_params:
            if not isinstance(strat_params[parameter], dict):
                strat_params[parameter] = {
                    'value': strat_params[parameter], 'vary': False}

            strata.add_par(parameter, **(strat_params[parameter]))

    return strata, sm, strata_type


def process_spatially_variable_strata(strata, component_name, yaml_data,
                                      locations, sm, strata_type):
    """
    This function handles the creation of stratigraphy components in
    openiam_cf.py in cases where the stratigraphy is spatially variable.

    Note that when unit depths decrease, the required changes in depth are
    accomplished by making the top unit (a shale) thicker. When unit depths
    decrease so much that the unit piches out (i.e., thickness in the
    subsurface goes to zero), the thickness will only go down to the minimum
    value of 1 m.

    :param component_name: name of the component model being handled. Note that
        the component name is used to extract the location data from locations
    :type component_name: str

    :param yaml_data: dictionary of input values from the .yaml file
    :type yaml_data: dict

    :param locations: dictionary linking specific components to the
        corresponding coordinates, number of locations, and component type
        (e.g., 'well': True)
    :type locations: dict

    :param sm: system model object
    :type sm: openiam.iam_base_classes.SystemModel

    :param strata_type: indicates whether the stratigraphy component is a
        'Stratigraphy', 'LookupTableStratigraphy', or 'DippingStratigraphy'
        component
    :type strata_type: str

    :returns: strata, sm
    """
    types_strata_obs = get_comp_types_strata_obs()

    comp_data = yaml_data[strata_type]

    # First, find the base name of the component (e.g., excluding '_001') and
    # the number within the name ('001'), which refers to the location.
    base_name = component_name[0:component_name.index('_')]
    comp_number = component_name[(component_name.index('_') + 1):None]
    comp_number = int(comp_number)

    if base_name in locations:
        # Get the coordinates for this component model
        comp_x_val = locations[base_name]['coordx'][comp_number]
        comp_y_val = locations[base_name]['coordy'][comp_number]

    elif 'DistributedTo' in yaml_data[base_name] and len(yaml_data[
            base_name]['DistributedTo']) != 0:
        # Get the coordinates from the component model connected to this one
        comp_x_val = locations[yaml_data[base_name]['DistributedTo'][0]]['coordx'][comp_number]
        comp_y_val = locations[yaml_data[base_name]['DistributedTo'][0]]['coordy'][comp_number]

    elif 'Connection' in yaml_data[base_name]:
        # Get the coordinates from the component model connected to this one
        comp_x_val = locations[yaml_data[base_name]['Connection']]['coordx'][comp_number]
        comp_y_val = locations[yaml_data[base_name]['Connection']]['coordy'][comp_number]

    if strata_type in types_strata_obs:
        comp_x_val = np.asarray(comp_x_val).flatten()
        comp_y_val = np.asarray(comp_y_val).flatten()

        if strata_type == 'LookupTableStratigraphy':
            # file_directory and file_name are handled in the connect_with_system() method
            strata.append(sm.add_component_model_object(LookupTableStratigraphy(
                name='strata' + component_name, parent=sm, intr_family='stratigraphy',
                locX=comp_x_val, locY=comp_y_val)))
        elif strata_type == 'DippingStratigraphy':
            strata.append(sm.add_component_model_object(DippingStratigraphy(
                name='strata' + component_name, parent=sm,
                locX=comp_x_val, locY=comp_y_val)))

            # This method sets the locXRef and locYRef from comp_data
            strata[-1].check_data_for_ref_loc(comp_data)

        if len(comp_x_val) == 1:
            strata[-1].grid_obs_keys = []

        if 'Parameters' in comp_data:
            strat_params = comp_data['Parameters']
        else:
            strat_params = comp_data

        for parameter in strat_params:
            if not parameter in ['ReferenceLocation', 'Controls', 'FileDirectory',
                                 'FileName']:
                if not isinstance(strat_params[parameter], dict):
                    strat_params[parameter] = {
                        'value': strat_params[parameter], 'vary': False}

                strata[-1].add_par(parameter, **(strat_params[parameter]))

    return strata, sm


def get_strata_type_from_yaml(yaml_data):
    """
    Function that obtains the type of stratigraphy component from a .yaml file.

    :param yaml_data: The dictionary created from an input .yaml_file
    :type yaml_data: dict

    :returns: comp_type
    """
    if 'Stratigraphy' in yaml_data:
        comp_type = 'Stratigraphy'

    elif 'LookupTableStratigraphy' in yaml_data:
        comp_type = 'LookupTableStratigraphy'

    elif 'DippingStratigraphy' in yaml_data:
        comp_type = 'DippingStratigraphy'

    else:
        err_msg = ''.join([
            'The code attempted to find the type of stratigraphy component being ',
            'used by checking the yaml_data dictionary. No applicable stratigraphy ',
            'component type was found, however. For the input dictionary to be ',
            'set up correctly, it needs to include one of the following keys: ',
            'Stratigraphy, DippingStratigraphy, or LookupTableStratigraphy. ',
            'Because no applicable component type was found, the code will proceed ',
            'with the Stratigraphy type. These conditions may lead to errors, ',
            'however. Check your input.'])
        logging.error(err_msg)

        comp_type = 'Stratigraphy'

    return comp_type


def get_default_num_shale_layers(strata_type, strata_data):
    """
    Although each stratigraphy component has a method that returns the number of 
    shale layers, that method is meant to be used after the component has been 
    set up (i.e., after the connect_with_system() method). This function returns 
    the number of shale layers for a component before the component has been set 
    up. The number is returned after checking the input provided for the stratigraphy 
    component in the yaml_data dictionary.
    """
    numberOfShaleLayers = DEFAULT_NUM_SHALE_LAYERS_DICT.get(
        strata_type, DEFAULT_NUM_SHALE_LAYERS)
    
    if 'Parameters' in strata_data:
        if 'numberOfShaleLayers' in strata_data['Parameters']:
            if isinstance(strata_data['Parameters']['numberOfShaleLayers'], dict):
                if 'value' in strata_data['Parameters']['numberOfShaleLayers']:
                    numberOfShaleLayers = strata_data[
                        'Parameters']['numberOfShaleLayers']['value']
            else:
                numberOfShaleLayers = strata_data['Parameters']['numberOfShaleLayers']
    
    return numberOfShaleLayers


def get_strata_info_from_component(strata_comp):
    """
    Function that obtains stratigraphy information from an already created
    stratigraphy component. The information obtained includes the number of
    shale layers, the thickness of each shale layer, the thickness of each
    aquifer, and the thickness of the reservoir. Minimum and maximum values for
    each value are also obtained, although those values are zero if the
    stratigraphy component does not have minimum and mazimum values.

    Note that this function deals with stratigraphy components that provide unit
    thicknesses and depths as parameters (e.g., Stratigraphy compomnent). This
    function is not appropriate for components that provide unit thicknesses and
    depths as observations (e.g., LookupTableStratigraphy).

    :param strata_comp: Stratigraphy component used to obtain information
    :type strata_comp: instance of Stratigraphy component class

    """
    # Initialize output dictionary
    strata_dict = dict()

    # Defaults for parameter numberOfShaleLayers
    strata_dict['numberOfShaleLayers_min'] = 0
    strata_dict['numberOfShaleLayers_max'] = 0
    strata_dict['numberOfShaleLayers_vary'] = False

    numShaleLayers = strata_comp.get_num_shale_layers()

    strata_dict['numberOfShaleLayers'] = numShaleLayers

    # The min, max, and numShaleLayers_vary are not used right now, but they
    # are kept here because they might be used in future versions.
    if 'numberOfShaleLayers' in strata_comp.pars:
        strata_dict['numberOfShaleLayers_min'] = \
            strata_comp.pars['numberOfShaleLayers'].min
        strata_dict['numberOfShaleLayers_max'] = \
            strata_comp.pars['numberOfShaleLayers'].max
        strata_dict['numberOfShaleLayers_vary'] = True

    # Initialize additional keys
    for key in ['shaleThicknesses', 'aquiferThicknesses',
                'shaleDepths', 'aquiferDepths',
                'shaleMidDepths', 'aquiferMidDepths',
                'shaleTopDepths', 'aquiferTopDepths']:
        strata_dict[key] = []

    for key in ['shaleThicknesses_min', 'shaleThicknesses_max',
                'aquiferThicknesses_min', 'aquiferThicknesses_max']:
        strata_dict[key] = numShaleLayers*[0]

    for key in ['shaleThicknesses_vary', 'aquiferThicknesses_vary']:
        strata_dict[key] = numShaleLayers*[False]

    # Get the unit thicknesses
    # Shale layers
    for shaleRef in range(numShaleLayers):
        shale_par_nm = 'shale{}Thickness'.format(shaleRef + 1)
        strata_dict['shaleThicknesses'].append(
            cfi_commons.get_parameter_val(strata_comp, shale_par_nm))

        # The min and max are not used currently, but they are kept for potential
        # updates in the future.
        if shale_par_nm in strata_comp.pars:
            strata_dict['shaleThicknesses_min'][shaleRef] = \
                strata_comp.pars[shale_par_nm].min
            strata_dict['shaleThicknesses_max'][shaleRef] = \
                strata_comp.pars[shale_par_nm].max
            strata_dict['shaleThicknesses_vary'][shaleRef] = True

        shale_par_nm = 'shale{}Depth'.format(shaleRef + 1)
        strata_dict['shaleDepths'].append(
            cfi_commons.get_parameter_val(strata_comp, shale_par_nm))

        shale_par_nm = 'shale{}MidDepth'.format(shaleRef + 1)
        strata_dict['shaleMidDepths'].append(
            cfi_commons.get_parameter_val(strata_comp, shale_par_nm))

        shale_par_nm = 'shale{}TopDepth'.format(shaleRef + 1)
        strata_dict['shaleTopDepths'].append(
            cfi_commons.get_parameter_val(strata_comp, shale_par_nm))

    # Aquifer layers
    for shaleRef in range(numShaleLayers-1):
        aq_par_nm = 'aquifer{}Thickness'.format(shaleRef + 1)

        strata_dict['aquiferThicknesses'].append(
            cfi_commons.get_parameter_val(strata_comp, aq_par_nm))

        if aq_par_nm in strata_comp.pars:
            strata_dict['aquiferThicknesses_min'][shaleRef] = \
                strata_comp.pars[aq_par_nm].min
            strata_dict['aquiferThicknesses_max'][shaleRef] = \
                strata_comp.pars[aq_par_nm].max
            strata_dict['aquiferThicknesses_vary'][shaleRef] = True

        aq_par_nm = 'aquifer{}Depth'.format(shaleRef + 1)
        strata_dict['aquiferDepths'].append(
            cfi_commons.get_parameter_val(strata_comp, aq_par_nm))

        aq_par_nm = 'aquifer{}MidDepth'.format(shaleRef + 1)
        strata_dict['aquiferMidDepths'].append(
            cfi_commons.get_parameter_val(strata_comp, aq_par_nm))

        aq_par_nm = 'aquifer{}TopDepth'.format(shaleRef + 1)
        strata_dict['aquiferTopDepths'].append(
            cfi_commons.get_parameter_val(strata_comp, aq_par_nm))

    # Add additional keys
    strata_dict['depth_default'] = False
    for par_name in ['reservoirThickness', 'depth', 'datumPressure']:
        for key in ['_min', '_max']:
            # Default values
            strata_dict[par_name+key] = 0  # e.g. 'reservoirThickness_min'

        strata_dict[par_name+'_vary'] = False   # e.g. 'reservoirThickness_vary'

        strata_dict[par_name] = cfi_commons.get_parameter_val(
            strata_comp, par_name)

        # The min and max vlaues are not used currently, but they are kept for
        # potential updates in the future.
        if par_name in strata_comp.pars:
            strata_dict[par_name+'_min'] = strata_comp.pars[par_name].min
            strata_dict[par_name+'_max'] = strata_comp.pars[par_name].max
            strata_dict[par_name+'_vary'] = True
        elif par_name in strata_comp.default_pars and par_name == 'depth':
            strata_dict['depth_default'] = True

    strata_dict['reservoirDepth'] = cfi_commons.get_parameter_val(
        strata_comp, 'reservoirDepth')

    strata_dict['reservoirMidDepth'] = cfi_commons.get_parameter_val(
        strata_comp, 'reservoirMidDepth')

    strata_dict['reservoirTopDepth'] = cfi_commons.get_parameter_val(
        strata_comp, 'reservoirTopDepth')

    return strata_dict


def get_unit_depth_from_component(numShaleLayers, stratigraphyComponent,
                                  unitNumber=1, unitType='shale',
                                  top_mid_bottom='top', depth_obs=False, sm=None):
    """
    This function  provides the top or bottom depth of a unit.

    Note that this function provides the unit depths that are derived from the
    simulation input, such as unit thicknesses, strike, and dip. In contrast,
    the function get_strata_info_from_component() deals only with the input provided
    directly to stratigraphy components (e.g., unit thicknesses and datumPressure).

    Unit depths can also be obtained with the composite parameters made during
    the connect_with_system() method of the Stratigraphy component. Those
    parameters work in forward simulations but are 0 in LHS and parstudy
    simulations, and this function solves that problem.

    :param numShaleLayers: Number of shale layers
    :type numShaleLayers: int

    :param stratigraphyComponent: Stratigraphy or LookupTableStratigraphy component
        used to obtain information
    :type stratigraphyComponent: instance of Stratigraphy or LookupTableStratigraphy
        component class

    :param unitNumber: The number for the unit targeted (e.g., 2 for aquifer 2)
    :type unitNumber: int

    :param unitType: The type of unit. This entry can only be 'shale', 'aquifer',
        or 'reservoir.'
    :type unitType: str

    :param top_mid_bottom: Option to obtain the 'top', 'mid', or 'bottom' depth
        of a unit
    :type top_mid_bottom: str

    :param depth_obs: Option specifying whether the type of stratigraphy component
        used produces unit depths as an observation. This is True for a
        LookupTableStratigraphy component.
    :type depth_obs: bool

    :param sm: system model object or None (default is None)
    :type sm: openiam.iam_base_classes.SystemModel or None

    :returns: unitDepth
    """
    if depth_obs:
        if unitType == 'shale' or unitType == 'aquifer':
            obs_nm = unitType + str(unitNumber)
        else:
            obs_nm = unitType

        if top_mid_bottom == 'bottom':
            obs_nm += 'Depth'
        elif top_mid_bottom == 'mid':
            obs_nm += 'MidDepth'
        elif top_mid_bottom == 'top':
            obs_nm += 'TopDepth'

        unitDepth = sm.collect_observations_as_time_series(
            stratigraphyComponent, obs_nm, indices=[0])[0]

    else:
        res_mid_check = False
        res_bottom_check = False

        # If the user wants the top depth of the reservoir, obtain the bottom of
        # shale 1 (deepest shale) instead.
        if unitType == 'reservoir':
            unitType = 'shale'
            unitNumber = 1

            # If the user wants the bottom depth of the reservoir, add the reservoir
            # thickness at the end.
            if top_mid_bottom == 'bottom':
                res_bottom_check = True
            elif top_mid_bottom == 'mid':
                res_mid_check = True
            elif top_mid_bottom == 'top':
                # The user wanted the top of the reservoir, which is the bottom of
                # shale 1
                top_mid_bottom = 'bottom'

        unitDepth = 0

        # First, handle the shallowest shale
        shale_par_nm = 'shale{}Thickness'.format(numShaleLayers)

        shaleThickness = cfi_commons.get_parameter_val(
            stratigraphyComponent, shale_par_nm)

        if unitNumber < numShaleLayers:
            unitDepth += shaleThickness  # start with thickness of the top shale
        else:  # unitNumber == numShaleLayers
            if unitType == 'shale' and top_mid_bottom == 'mid':
                unitDepth += (shaleThickness / 2)
            elif unitType == 'shale' and top_mid_bottom == 'bottom':
                unitDepth += shaleThickness

        if unitNumber != numShaleLayers:
            for shaleRef in range(numShaleLayers - 1, unitNumber - 1, -1):
                aq_par_nm = 'aquifer{}Thickness'.format(shaleRef)
                shale_par_nm = 'shale{}Thickness'.format(shaleRef)

                aquiferThickness = cfi_commons.get_parameter_val(
                    stratigraphyComponent, aq_par_nm)

                shaleThickness = cfi_commons.get_parameter_val(
                    stratigraphyComponent, shale_par_nm)

                if shaleRef > unitNumber:
                    unitDepth += (aquiferThickness + shaleThickness)

                elif shaleRef == unitNumber and unitType == 'shale':
                    unitDepth += aquiferThickness
                    if top_mid_bottom == 'mid':
                        unitDepth += (shaleThickness / 2)
                    elif top_mid_bottom == 'bottom':
                        unitDepth += shaleThickness
                elif shaleRef == unitNumber and unitType == 'aquifer':
                    if top_mid_bottom == 'mid':
                        unitDepth += (aquiferThickness / 2)
                    elif top_mid_bottom == 'bottom':
                        unitDepth += aquiferThickness

        if res_bottom_check or res_mid_check:
            par_name = 'reservoirThickness'

            reservoirThickness = cfi_commons.get_parameter_val(
                stratigraphyComponent, par_name)

            if res_mid_check:
                unitDepth +=( reservoirThickness / 2)
            elif res_bottom_check:
                unitDepth += reservoirThickness

    if not isinstance(unitDepth, float):
        unitDepth = float(unitDepth)

    return unitDepth


def get_strat_param_dict_for_link(sparam, strat_comp):
    """
    Checks for a given parameter name (sparam) in the composite, stochastic,
    deterministic, and default parameters of a stratigraphy type component. If
    the parameter name is found in one type of parameter, that parameter
    dictionary type is returned.
    """
    connect = None

    if sparam in strat_comp.composite_pars:
        connect = strat_comp.composite_pars
    elif sparam in strat_comp.pars:
        connect = strat_comp.pars
    elif sparam in strat_comp.deterministic_pars:
        connect = strat_comp.deterministic_pars
    elif sparam in strat_comp.default_pars:
        connect = strat_comp.default_pars
    else:
        info_msg = 'Unable to find parameter {} for component {}.'.format(
            sparam, strat_comp.name)
        logging.info(info_msg)

    return connect


def check_color_alpha_label_yaml_input(yaml_input, strata_plot_data, name):
    """
    Function that checks if color, alpha, or labels provided for stratigraphic
    units in a .yaml file are acceptable. If the input is recognized and deemed
    acceptable, the input is added to yaml_input. If not, the input is excluded
    from yaml_input and a warning message is logged and printed.

    :param yaml_input: Dictionary of keys provided for the plotting function
        used (e.g., stratigraphy_plot() or stratigraphic_column).
    :type yaml_input: dict

    :param strata_plot_data: Input values provided for the plot in the Plots
        section of a .yaml file.
    :type strata_plot_data: dict

    :param name: Name of the plot
    :type name: str
    """
    alphaNames = ['ReservoirAlpha', 'ShaleAlpha', 'AquiferAlpha', 'WellAlpha']

    alphaNames += ['Shale' + str(num) + 'Alpha' for num in range(1, 31)]
    alphaNames += ['Aquifer' + str(num) + 'Alpha' for num in range(1, 30)]

    colorNames = ['ReservoirColor', 'ShaleColor', 'AquiferColor', 'WellColor']

    colorNames += ['Shale' + str(num) + 'Color' for num in range(1, 31)]
    colorNames += ['Aquifer' + str(num) + 'Color' for num in range(1, 30)]

    labelNames = ['ReservoirLabel', 'WellLabel']

    labelNames += ['Shale' + str(num) + 'Label' for num in range(1, 31)]
    labelNames += ['Aquifer' + str(num) + 'Label' for num in range(1, 30)]

    for input_type in alphaNames:
        if input_type in strata_plot_data:
            try:
                if 0 < strata_plot_data[input_type] <= 1:
                    # The alpha input corresponds with the alpha used for filled areas
                    yaml_input[input_type + 'Fill'] = strata_plot_data[input_type]

                    # The alpha value used for lines is higher
                    yaml_input[input_type] = strata_plot_data[input_type] + (
                        (1 - strata_plot_data[input_type]) * ALPHA_SCALE_FACTOR)
                else:
                    color_alpha_label_error_msg(
                        input_type, strata_plot_data[input_type], 'alpha', name)
            except:
                color_alpha_label_error_msg(
                    input_type, strata_plot_data[input_type], 'alpha', name)

    for input_type in colorNames:
        if input_type in strata_plot_data:
            if isinstance(strata_plot_data[input_type], str):
                if is_color_like(strata_plot_data[input_type]):
                    yaml_input[input_type] = strata_plot_data[input_type]
                else:
                    color_alpha_label_error_msg(
                        input_type, strata_plot_data[input_type], 'str', name)

            elif isinstance(strata_plot_data[input_type], list):
                if not is_color_like(strata_plot_data[input_type]):
                    color_alpha_label_error_msg(
                        input_type, strata_plot_data[input_type], 'list', name)
                else:
                    yaml_input[input_type] = strata_plot_data[input_type]
            else:
                color_alpha_label_error_msg(
                    input_type, strata_plot_data[input_type], 'neither', name)

    for input_type in labelNames:
        if input_type in strata_plot_data:
            if isinstance(strata_plot_data[input_type], str):
                yaml_input[input_type] = strata_plot_data[input_type]
            else:
                color_alpha_label_error_msg(
                    input_type, strata_plot_data[input_type], 'label', name)

    return yaml_input


def color_alpha_label_error_msg(input_type, input_value, error_type, name):
    """
    Function that logs and prints an error message when invalid input is
    provided for a unit color, unit alpha, or unit label.
    """
    err_msg_pt1 = ''.join(['The ', input_type, ' input provided for the figure ',
                           name, '(', input_type, ':', ' ', input_value, ')'])
    if error_type == 'alpha':
        err_msg_pt2 = ''.join([' was not a value between 0 and 1.'])
    elif error_type == 'str':
        err_msg_pt2 = ''.join([' was a string input, but not one that matplotlib ',
                               'recognized as a color (e.g., "r", "teal", ',
                               'or "#030764").'])
    elif error_type == 'list':
        err_msg_pt2 = ''.join([' was a list, but not a list of length three ',
                               'that matplotlib recognized as a color (e.g., ',
                               '"[1, 0, 0]" or "[0.67, 0.33, 0]"). The three ',
                               'values in the list must be between 0 and 1 ',
                               '(for fractions of [Red, Green, Blue]).'])
    elif error_type == 'neither':
        err_msg_pt2 = ''.join([' was neither a string nor a list that matplotlib ',
                               'recognized as a color (e.g., "r", "teal", "#030764", ',
                               '"[1, 0, 0]," or "[0.67, 0.33, 0]"). If given as ',
                               'a list, the three values in the list must be ',
                               'between 0 and 1 (for fractions of [Red, Green, Blue]).'])
    elif error_type == 'label':
        err_msg_pt2 = ''.join([' was not a string. Only string inputs (e.g., ',
                               '"AZMI", "BGWP", or "Freshwater Aquifer") can be ',
                               'provided for unit labels. '])

    err_msg_pt3 = ' The figure {} will use the default approach for setting {}.'.format(
        name, input_type)

    err_msg = err_msg_pt1 + err_msg_pt2 + err_msg_pt3

    logging.debug(err_msg)


def get_plotting_setup_for_units(strat_plot_yaml_input, numShaleLayers,
                                 reservoirThickness=None, shaleThicknessList=None,
                                 aquiferThicknessList=None, include_thickness=False):
    """
    Function that provides the colors, alphas, and labels for each unit. The
    input provided in strat_plot_yaml_input is used if present. Otherwise,
    default values are used.

    :param strata_plot_data: Dictionary containing input values for a plot.
        This dictionary should already have been updated with the function
        check_color_alpha_label_yaml_input().
    :type strata_plot_data: dict

    :param numShaleLayers: Number of shale layers.
    :type numShaleLayers: int

    :param reservoirThickness: Thickness of the reservoir (m)
    :type reservoirThickness: int or float

    :param shaleThicknessList: List of shale thicknesses (m)
    :type shaleThicknessList: list

    :param aquiferThicknessList: List of aquifer thicknesses (m)
    :type aquiferThicknessList: list
    """
    reservoirColor = strat_plot_yaml_input.get(
        'ReservoirColor', UNIT_COLOR_DICT['ReservoirColor'])
    reservoirAlpha = strat_plot_yaml_input.get(
        'ReservoirAlpha', UNIT_ALPHA_DICT['ReservoirAlpha'])
    reservoirAlphaFill = strat_plot_yaml_input.get(
        'ReservoirAlphaFill', UNIT_ALPHA_FILL_DICT['ReservoirAlphaFill'])
    if include_thickness and not reservoirThickness is None:
        reservoirLabel = strat_plot_yaml_input.get(
            'ReservoirLabel', UNIT_LABEL_THICKNESS_DICT['ReservoirLabel'].format(
                reservoirThickness))
    else:
        reservoirLabel = strat_plot_yaml_input.get(
            'ReservoirLabel', UNIT_LABEL_DICT['ReservoirLabel'])

    wellColor = strat_plot_yaml_input.get('WellColor', WELL_COLOR)
    wellAlpha = strat_plot_yaml_input.get('WellAlpha', WELL_ALPHA)
    wellAlphaFill = strat_plot_yaml_input.get('WellAlphaFill', WELL_ALPHA_FILL)
    wellLabel = strat_plot_yaml_input.get('WellLabel', None)

    shaleColor = [None] * numShaleLayers
    shaleAlpha = [None] * numShaleLayers
    shaleAlphaFill = [None] * numShaleLayers
    shaleLabel = [None] * numShaleLayers

    aquiferColor = [None] * (numShaleLayers - 1)
    aquiferAlpha = [None] * (numShaleLayers - 1)
    aquiferAlphaFill = [None] * (numShaleLayers - 1)
    aquiferLabel = [None] * (numShaleLayers - 1)

    for shaleRef in range(0, numShaleLayers):
        shaleColor[shaleRef] = strat_plot_yaml_input.get(
            'ShaleColor', UNIT_COLOR_DICT['ShaleColor'])

        # More specific color input overrides the more general color input
        if 'Shale{:.0f}Color'.format(shaleRef + 1) in strat_plot_yaml_input:
            shaleColor[shaleRef] = strat_plot_yaml_input[
                'Shale{:.0f}Color'.format(shaleRef + 1)]

        shaleAlpha[shaleRef] = strat_plot_yaml_input.get(
            'ShaleAlpha', UNIT_ALPHA_DICT['ShaleAlpha'])

        shaleAlphaFill[shaleRef] = strat_plot_yaml_input.get(
            'ShaleAlphaFill', UNIT_ALPHA_FILL_DICT['ShaleAlphaFill'])

        # More specific alpha input overrides the more general color input
        if 'Shale{:.0f}Alpha'.format(shaleRef + 1) in strat_plot_yaml_input:
            shaleAlpha[shaleRef] = strat_plot_yaml_input[
                'Shale{:.0f}Alpha'.format(shaleRef + 1)]

        if 'Shale{:.0f}AlphaFill'.format(shaleRef + 1) in strat_plot_yaml_input:
            shaleAlphaFill[shaleRef] = strat_plot_yaml_input[
                'Shale{:.0f}AlphaFill'.format(shaleRef + 1)]

        if include_thickness and not shaleThicknessList is None:
            shaleLabel[shaleRef] = strat_plot_yaml_input.get(
                'Shale{:.0f}Label'.format(shaleRef + 1),
                UNIT_LABEL_THICKNESS_DICT['ShaleLabel'].format(
                    shaleRef + 1, shaleThicknessList[shaleRef]))
        else:
            shaleLabel[shaleRef] = strat_plot_yaml_input.get(
                'Shale{:.0f}Label'.format(shaleRef + 1),
                UNIT_LABEL_DICT['ShaleLabel'].format(shaleRef + 1))

        # Get input for aquifers
        if shaleRef != (numShaleLayers - 1):
            aquiferColor[shaleRef] = strat_plot_yaml_input.get(
                'Aquifer{:.0f}Color'.format(shaleRef + 1), UNIT_COLOR_DICT['AquiferColor'])

            # More specific color input overrides the more general color input
            if 'Aquifer{:.0f}Color'.format(shaleRef + 1) in strat_plot_yaml_input:
                aquiferColor[shaleRef] = strat_plot_yaml_input[
                    'Aquifer{:.0f}Color'.format(shaleRef + 1)]

            aquiferAlpha[shaleRef] = strat_plot_yaml_input.get(
                'AquiferAlpha', UNIT_ALPHA_DICT['AquiferAlpha'])

            aquiferAlphaFill[shaleRef] = strat_plot_yaml_input.get(
                'AquiferAlphaFill', UNIT_ALPHA_FILL_DICT['AquiferAlphaFill'])

            # More specific alpha input overrides the more general color input
            if 'Aquifer{:.0f}Alpha'.format(shaleRef + 1) in strat_plot_yaml_input:
                aquiferAlpha[shaleRef] = strat_plot_yaml_input[
                    'Aquifer{:.0f}Alpha'.format(shaleRef + 1)]

            if 'Aquifer{:.0f}AlphaFill'.format(shaleRef + 1) in strat_plot_yaml_input:
                aquiferAlphaFill[shaleRef] = strat_plot_yaml_input[
                    'Aquifer{:.0f}AlphaFill'.format(shaleRef + 1)]

            if include_thickness and not aquiferThicknessList is None:
                aquiferLabel[shaleRef] = strat_plot_yaml_input.get(
                    'Aquifer{:.0f}Label'.format(shaleRef + 1),
                    UNIT_LABEL_THICKNESS_DICT['AquiferLabel'].format(
                        shaleRef + 1, aquiferThicknessList[shaleRef]))
            else:
                aquiferLabel[shaleRef] = strat_plot_yaml_input.get(
                    'Aquifer{:.0f}Label'.format(shaleRef + 1),
                    UNIT_LABEL_DICT['AquiferLabel'].format(shaleRef + 1))

    return reservoirColor, reservoirAlpha, reservoirAlphaFill, reservoirLabel, \
        shaleColor, shaleAlpha, shaleAlphaFill, shaleLabel, \
            aquiferColor, aquiferAlpha, aquiferAlphaFill, aquiferLabel, \
                wellColor, wellAlpha, wellAlphaFill, wellLabel
