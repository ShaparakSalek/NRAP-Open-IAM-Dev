import os
import sys
import re
import logging
import numpy as np

try:
    from openiam.components.iam_base_classes import IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cf_interface.strata import (get_strata_info_from_component,
                                         get_unit_depth_from_component,
                                         get_comp_types_strata_pars,
                                         get_comp_types_strata_obs)


def process_parameters(component, component_data, name2obj_dict=None):
    """
    Check whether any parameters are provided and process their data.

    Process parameters of the component based on the information
    provided in the component data.

    One can provide input to a parameter in the .yaml file as a string.
    These strings must correspond to the bottom depth of a unit, such as
    shale1Depth or aquifer2Depth. For example, one might want to set the wellTop
    or reservoirDepth parameters for an OpenWellbore in this way. Note that
    reservoirDepth can be represented by shale1Depth. Using these string inputs
    is especially useful when the strike and dip option is used for the
    stratigraphy, such that there are spatial variations in unit depths. The
    code obtains the stratigraphy component created for each component, so the
    depths can vary as needed.
    """
    types_strata_pars = get_comp_types_strata_pars()
    types_strata_obs = get_comp_types_strata_obs()

    if ('Parameters' in component_data) and (component_data['Parameters']):
        for key in component_data['Parameters']:
            
            # The case for "par_name: par_value"
            if not isinstance(component_data['Parameters'][key], dict):
                component_data['Parameters'][key] = {
                    'value': component_data['Parameters'][key], 'vary': False}

            if name2obj_dict is not None and 'value' in component_data['Parameters'][key]:
                value = component_data['Parameters'][key]['value']

                # Get the stratigraphy component
                strata_comp = name2obj_dict['strata']

                if isinstance(value, str) and name2obj_dict['strata_type'] in types_strata_pars:
                    strata_dict = get_strata_info_from_component(strata_comp)

                    numShaleLayers = strata_dict['numberOfShaleLayers']
                    shaleThicknesses = strata_dict['shaleThicknesses']
                    aquiferThicknesses = strata_dict['aquiferThicknesses']

                    unitType, unitNum, metricType = handle_str_input(
                        numShaleLayers, value, component.name, key)

                    if metricType == 'Thickness':
                        if unitType == 'shale':
                            thicknessValue = shaleThicknesses[unitNum - 1]

                        elif unitType == 'aquifer':
                            thicknessValue = aquiferThicknesses[unitNum - 1]

                        component_data['Parameters'][key]['value'] = thicknessValue

                    elif metricType == 'Depth':
                        if 'TopDepth' in value:
                            top_mid_bottom = 'top'
                        elif 'MidDepth' in value:
                            top_mid_bottom = 'mid'
                        else:
                            top_mid_bottom = 'bottom'

                        depthValue = get_unit_depth_from_component(
                            numShaleLayers, strata_comp, unitNumber=unitNum,
                            unitType=unitType, top_mid_bottom=top_mid_bottom)

                        component_data['Parameters'][key]['value'] = depthValue
                    
                    check_min_max_vals(component_data['Parameters'], key, component.name)

                    component.add_par(key, **component_data['Parameters'][key])

                    # Reset the value so that other locations for the same type
                    # of component will be handled properly.
                    component_data['Parameters'][key]['value'] = value

                elif isinstance(value, str) and name2obj_dict[
                        'strata_type'] in types_strata_obs:
                    component.add_par_linked_to_obs(key, strata_comp.linkobs[value])

                else:
                    check_min_max_vals(component_data['Parameters'], key, component.name)
                    component.add_par(key, **component_data['Parameters'][key])

            else:
                check_min_max_vals(component_data['Parameters'], key, component.name)
                
                component.add_par(key, **component_data['Parameters'][key])


def check_min_max_vals(param_dict, key, comp_name):
    """
    Checks if the parameter is stochastic. If so, checks if the minimum value is 
    greater than or equal to the maximum value. If so, raises an error.
    """
    stochastic_check = ('min' in param_dict[key] and 'max' in param_dict[key])
    
    if stochastic_check:
        min_val = param_dict[key]['min']
        max_val = param_dict[key]['max']
        
        if min_val >= max_val:
            err_msg = ''.join([
                'For the parameter {} of the component {}, the minimum and ', 
                'maximum values were given as {} and {}, respectively. The ', 
                'minimum value cannot be greater than or equal to the maximum ', 
                'value. Check your input.']).format(
                    key, comp_name, min_val, max_val)
            
            logging.error(err_msg)
            raise KeyError(err_msg)


def process_dynamic_inputs(component, component_data, array_like=False,
                           check_second_dim=False, **kwargs):
    """
    Check whether any dynamic inputs are provided and process their data.

    Process dynamic inputs of the component based on the information
    provided in the component data.

    component - object of ComponentModel class
    component_data - dictionary containing data about setup of the component
    array_like - flag showing whether slice of the data provided through
    dynamic input setup of the component should be array like

    if array_like is True and check_second_dim is True as well, it allows to check whether
    the slice of the data (second dimension of the data) has the right size.

    In kwargs:
        size
        quantity_to_compare
        For size, quantity_to_compare see documentation for
        check_size method
    """
    if ('DynamicParameters' in component_data) and (component_data['DynamicParameters']):
        for key in component_data['DynamicParameters']:
            if key != 'structure':
                inp_data = component_data['DynamicParameters'][key]
                if isinstance(inp_data, str):  # if filename is provided
                    inp_data_file_path = os.path.join(IAM_DIR, inp_data)
                    if os.path.isfile(inp_data_file_path):
                        data = np.genfromtxt(inp_data_file_path,
                                             delimiter=",", dtype='f8',
                                             comments='#')
                        if data.ndim == 2:
                            data = data.T
                        else:
                            if array_like:
                                data = data.reshape(data.shape[0], -1)

                    else:
                        err_msg = '{} is not a valid file name.'.format(inp_data)
                        logging.error(err_msg)
                        raise FileNotFoundError(err_msg)
                else:  # if numerical data is provided
                    if array_like:
                        data = np.array(inp_data)
                        data = data.reshape(data.shape[0], -1)
                    else:
                        data = inp_data

                if check_second_dim:
                    check_size(component, data, kwargs['size'],
                               kwargs['quantity_to_compare'], key)
                component.add_dynamic_kwarg(key, data)


def check_size(component, data, size, quantity_to_compare, dyn_input_name):
    """
    data is a matrix like np.array whose size needs to be checked
    size is a required size of a slice
    quantity_to_compare is a string indicating name of a quantity (if applicable)
    whose size should be equal to the second dimension of the data
    dyn_input_name is a string indicating name of a dynamic input (usually,
    the one used in the model method of component, e.g.,
    pressure/CO2saturation/co2rate, etc.)
    """
    if data.shape[1] != size:
        err_msg = ''.join([
            '{} does not coincide with the second dimension of the provided ',
            'dynamic keyword argument {} for component {}']).format(
                quantity_to_compare, dyn_input_name, component.name)
        logging.error(err_msg)
        raise IndexError(err_msg)


def handle_str_input(num_shale_layers, str_input, cmpnt_name, par_name):
    """
    This function takes string input provided for a parameter in a .yaml file
    and returns the corresponding unitType ('shale' or 'aquifer'), unitNum
    (unit number, 1, 2, 3, etc.), and metricType ('Thickness', 'Depth', 'MidDepth',
    or 'TopDepth').

    :param num_shale_layers: number of shale layers
    :type num_shale_layers: int

    :param str_input: input for parameter provided instead of value
    :type str_input: str

    :param cmpnt_name: name of the component object for which the given parameter
        is processed
    :type cmpnt_name: str

    :param par_name: name of the parameter for which the string input is
        provided
    :type par_name: str
    """
    # Create a list of possible stratigraphy parameters or observations
    potential_str_references = [
        'shale{}Depth'.format(ind) for ind in range(1, num_shale_layers+1)]+[
            'shale{}MidDepth'.format(ind) for ind in range(1, num_shale_layers+1)]+[
                'shale{}TopDepth'.format(ind) for ind in range(1, num_shale_layers+1)]+[
                    'shale{}Thickness'.format(ind) for ind in range(1, num_shale_layers+1)]

    potential_str_references += [
        'aquifer{}Depth'.format(ind) for ind in range(1, num_shale_layers+1)]+[
            'aquifer{}MidDepth'.format(ind) for ind in range(1, num_shale_layers+1)]+[
                'aquifer{}TopDepth'.format(ind) for ind in range(1, num_shale_layers+1)]+[
                    'aquifer{}Thickness'.format(ind) for ind in range(1, num_shale_layers+1)]

    potential_str_references += [
        'reservoirDepth', 'reservoirMidDepth', 'reservoirTopDepth', 'reservoirThickness']

    # If string provided as parameter value is among possible stratigraphy parameters
    # determine whether it's for shale or aquifer and whether it's depth or thickness
    if str_input in potential_str_references:
        if 'reservoir' in str_input:
            unitType = 'reservoir'
            unitNum = 0
            metricType = str_input[9:]
        else:
            # Match the str_input to the pattern
            out = re.match(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", str_input)
            unitType = out.group(1)
            unitNum = int(out.group(2))
            metricType = out.group(3)

    else:
        err_msg = ''.join([
            'A string ({}) was used as input for the parameter {} in the component {}. ',
            'The string did not match any of the expected values, however. ',
            'Expected values are the depths to the bottom, middle, or top of a unit ',
            '(e.g., shale1Depth, aquifer2MidDepth, shale3TopDepth, or reservoirDepth). ',
            'Check the input in the .yaml file.']).format(str_input, par_name, cmpnt_name)
        raise KeyError(err_msg)

    return unitType, unitNum, metricType


def get_parameter_val(comp, par_name, sm=None, yaml_data=None):
    """
    Function that checks if a parameter (par_name) is in the stochastic,
    deterministic, composite, or default parameters of a component. If so, the
    parameter value is returned.
    """
    par_val = None

    # If the parameter is a composite parameter that returns as zero (which
    # can happen for lhs or parstudy simulations) or linked to another parameter,
    # the function check_parameter_types() updates run_again_keys and then runs again.
    # When it runs again, it is directed to another component from which to obtain
    # the parameter value (e.g., getting 'aquifer2Depth' from a stratigraphy
    # component in lieu of 'wellTop').
    run_again_keys = [
        'comp', 'par_name', 'unit_number', 'unit_type', 'top_mid_bottom']

    run_again_info = {key: None for key in run_again_keys}

    par_val, run_again, run_again_info = check_parameter_types(
        comp, par_name, run_again_info, sm=sm, yaml_data=yaml_data)

    if run_again:
        par_val, run_again, run_again_info = check_parameter_types(
            run_again_info['comp'], run_again_info['par_name'], run_again_info,
            sm=sm, yaml_data=yaml_data)

    return par_val


def check_parameter_types(comp, par_name, run_again_info, sm=None, yaml_data=None):
    """
    This function is used by get_parameter_val() to check if the parameter in
    question (par_name) is in the stochastic (.pars), deterministic, composite,
    linked, or default parameters. It will return the value once the parameter
    is found. If a composite parameter returns as 0 (which happens in lhs and parstudy
    simulations) or if the parameter is linked to another parameter, this function
    provides the component and parameter name that will provide the correct value
    (e.g., 'aquifer2Depth' for the 'wellTop' parameter).
    """
    # These lists indicate the stratigraphy component types that offer thicknesses
    # and depths as parameters or as observations.
    types_strata_pars = get_comp_types_strata_pars()
    types_strata_obs = get_comp_types_strata_obs()

    run_again = False

    par_val = None

    # This list is only meant to include composite parameters, not parameters
    # linked to other parameters (e.g., aqu_thick or top_depth for a GenericAquifer).
    strat_related_pars = ['wellTop', 'reservoirDepth', 'wellDepth']

    if par_name in comp.pars:
        par_val = comp.pars[par_name].value

    elif par_name in comp.deterministic_pars:
        par_val = comp.deterministic_pars[par_name].value

    elif par_name in comp.obslinked_pars:
        par_connection = comp.obslinked_pars[par_name]

        linked_comp_name = par_connection[0:par_connection.index('.')]

        run_again_info['comp'] = sm.component_models[linked_comp_name]

        run_again = True

        if par_name == 'wellTop':
            LeakTo = yaml_data[comp.name]['LeakTo']
            run_again_info['unit_number'] = int(LeakTo[7:])
            run_again_info['top_mid_bottom'] = 'bottom'
            run_again_info['par_name'] = LeakTo + 'Depth'

        elif par_name in ['reservoirDepth', 'wellDepth']:
            run_again_info['par_name'] = 'shale1Depth'
            run_again_info['unit_number'] = 1
            run_again_info['top_mid_bottom'] = 'bottom'

    elif 'Depth' in par_name and (comp.class_type in types_strata_obs):
        run_again_info = get_strat_run_again_info(
            par_name, run_again_info)

        par_val = get_unit_depth_from_component(
            comp.get_num_shale_layers(), comp,
            unitNumber=run_again_info['unit_number'],
            unitType=run_again_info['unit_type'],
            top_mid_bottom=run_again_info['top_mid_bottom'],
            depth_obs=True, sm=sm)

    elif par_name in comp.composite_pars:
        # In cases where a parameter is set as a composite parameter (e.g.,
        # setting the wellTop or reservoirDepth parameter of an OpenWellore
        # only by having 'LeakTo: aquifer2') and the simulation uses lhs or
        # parstudy analysis types, the composite_pars value will be returned
        # as zero. The composite parameters are handled correctly during the
        # simulation, but attempting to access the value afterwards (i.e., while
        # processing the plots) returns a zero value. The 'if' statements below
        # address such situations - additional cases may need to be added for
        # new component parameters.
        par_val = comp.composite_pars[par_name].value

        if par_val == 0:
            if 'Depth' in par_name and (comp.class_type in types_strata_pars):
                run_again_info = get_strat_run_again_info(
                    par_name, run_again_info)

                par_val = get_unit_depth_from_component(
                    comp.get_num_shale_layers(), comp,
                    unitNumber=run_again_info['unit_number'],
                    unitType=run_again_info['unit_type'],
                    top_mid_bottom=run_again_info['top_mid_bottom'])

            elif par_name in strat_related_pars and not yaml_data is None:

                try:
                    run_again_info['comp'] = sm.component_models['strata']
                except KeyError:
                    # Case for spatially variable stratigraphy
                    run_again_info['comp'] = sm.component_models['strata' + comp.name]

                if par_name == 'wellTop':
                    LeakTo = yaml_data[comp.name]['LeakTo']
                    run_again_info['unit_number'] = int(LeakTo[7:])
                    run_again_info['top_mid_bottom'] = 'bottom'
                    run_again_info['par_name'] = LeakTo + 'Depth'

                elif par_name in ['reservoirDepth', 'wellDepth']:
                    run_again_info['par_name'] = 'shale1Depth'
                    run_again_info['unit_number'] = 1
                    run_again_info['top_mid_bottom'] = 'bottom'

                run_again = True

            else:
                # This is where cases for new stratigraphy-related parameters should be entered
                pass

    elif par_name in comp.parlinked_pars and sm is not None:
        par_connection = comp.parlinked_pars[par_name]

        linked_comp_name = par_connection[0:par_connection.index('.')]
        run_again_info['par_name'] = par_connection[(par_connection.index('.') + 1):None]

        run_again_info['comp'] = sm.component_models[linked_comp_name]

        run_again = True

    else:
        err_msg = ''.join([
            'The parameter {} was not found for the component {}. This issue ',
            'can occur, for example, if the parameter was not added as a composite ',
            'parameter of the {} component being used. Check your input.'
            ]).format(par_name, comp.name, comp.class_type)

        if 'Thickness' in par_name and ('shale' in par_name or 'aquifer' in par_name):
            # Parameters like 'shale1Thickness' will not be in the default parameters,
            # but 'shaleThickness' will be.
            if 'shale' in par_name:
                par_name = 'shaleThickness'
            elif 'aquifer' in par_name:
                par_name = 'aquiferThickness'
        elif 'Depth' in par_name:
            logging.error(err_msg)

        try:
            par_val = comp.default_pars[par_name].value
        except:
            raise KeyError(err_msg)

    return par_val, run_again, run_again_info


def get_strat_run_again_info(par_name, run_again_info):
    """
    Updates the run_again_info dictionary for check_parameter_types() in cases
    where a value is needed from a stratigraphy component (e.g., Stratigraphy,
    LookupTableStratigraphy, or DippingStratigraphy). The run_again_info
    dictionary is then used to obtain a value with the get_unit_depth_from_component()
    function.
    """
    if 'aquifer' in par_name:
        run_again_info['unit_type'] = 'aquifer'
        par_name_end = par_name[7:]
    elif 'shale' in par_name:
        run_again_info['unit_type'] = 'shale'
        par_name_end = par_name[5:]
    elif 'reservoir' in par_name:
        run_again_info['unit_type'] = 'reservoir'
        par_name_end = par_name[9:]
        run_again_info['unit_number'] = 0

    if 'TopDepth' in par_name_end or 'MidDepth' in par_name_end:
        if run_again_info['unit_type'] != 'reservoir':
            run_again_info['unit_number'] = int(par_name_end[:-8])

        if 'TopDepth' in par_name_end:
            run_again_info['top_mid_bottom'] = 'top'
        elif 'MidDepth' in par_name_end:
            run_again_info['top_mid_bottom'] = 'mid'

    elif 'Depth' in par_name_end:
        if run_again_info['unit_type'] != 'reservoir':
            run_again_info['unit_number'] = int(par_name_end[:-5])

        run_again_info['top_mid_bottom'] = 'bottom'

    return run_again_info
