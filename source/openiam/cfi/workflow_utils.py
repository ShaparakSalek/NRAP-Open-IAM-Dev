
import logging
import openiam as iam
import numpy as np

YAML_INPUT_WARNING_MSG = ''.join([
    'The input provided for {} under the Workflow section of the .yaml file ',
    'was not of type {}. The default setting of {} will be used.'])


def get_automation_input(yaml_data):
    """
    Checks for input related to the automatic setup of a .yaml control file for
    a Workflow.
    """
    automation_input = {}
    for key in ['AutomateResCompSetup', 'AutomateWellCompSetup',
                'AutomateAqCompSetup', 'AutomatePlotsSetup']:
        inp_key = key[8:-5]
        automation_input[inp_key] = True
        if key in yaml_data['Workflow']['Options']:
            if isinstance(yaml_data['Workflow']['Options'][key], bool):
                automation_input[inp_key] = yaml_data['Workflow']['Options'][key]
            else:
                warning_msg = YAML_INPUT_WARNING_MSG.format(key, 'boolean', True)
                logging.warning(warning_msg)

    return automation_input


def get_parameter_names(component_type):
    """
    Function that returns the names of all parameters for a given component type.
    """
    num_years = 5
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    dummy_sm = iam.SystemModel(model_kwargs=sm_model_kwargs)

    comp_type = getattr(iam, component_type)

    dummy_comp = dummy_sm.add_component_model_object(comp_type(
        name='dummy_comp', parent=dummy_sm))

    par_names = dummy_comp.default_pars.keys()

    return par_names
