import os
import sys
import logging

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)

try:
    from openiam import IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

# import workflow utilities
from openiam.cfi.workflow_utils import get_parameter_names, get_automation_input

# Can be expanded to incorporate new workflows
from openiam.cfi.workflows import AOR_workflow

# Can be expanded to incorporate new workflows
WORKFLOW_OPTIONS = ['AoR']

YAML_INPUT_WARNING_MSG = ''.join([
    'The input provided for {} under the Workflow section of the .yaml file ',
    'was not of type {}. The default setting of {} will be used.'])


def iam_workflow_setup(yaml_data, strata):
    """
    This function reads the input provided in the Workflow section of a .yaml
    file. Depending on the options provided, the function can add to the yaml_data
    dictionary which is then returned by the function. For example, The AoR
    Workflow requires the inclusion of AoR plot entries using multiple metrics
    (e.g., pressure, CO2saturation, Dissolved_salt_volume, and Dissolved_CO2_volume).
    Setting up these plots and the required components can require an experienced
    user, so this function is designed to handle much of that effort.
    """

    if 'Type' in yaml_data['Workflow']:
        if yaml_data['Workflow']['Type'] in WORKFLOW_OPTIONS:
            workflow_type = yaml_data['Workflow']['Type']
        else:
            warning_msg = ''.join([
                'A Workflow section was included in the .yaml file, but the Workflow ',
                'Type provided (', yaml_data['Workflow']['Type'], ') was not ',
                'recognized as one of the available options (', WORKFLOW_OPTIONS,
                '). Therefore, the input in the Workflow section will not be used.'])
            logging.warning(warning_msg)
    else:
        info_msg = ''.join([
            'A Workflow section was included in the .yaml file, but the section ',
            'did not include a Type entry. Therefore, the input in the Workflow ',
            'section will not be used.'])
        logging.info(info_msg)

    if 'Options' not in yaml_data['Workflow']:
        info_msg = ''.join([
            'The Workflow section of the .yaml file did not include an ',
            'Options section. The default settings will be used.'])
        logging.info(info_msg)

        yaml_data['Workflow']['Options'] = {}


    if workflow_type == 'AoR':
        automated_options = ['AutomateResCompSetup', 'AutomateWellCompSetup',
                             'AutomateAqCompSetup', 'AutomatePlotsSetup']
        automation_input = get_automation_input(yaml_data, automated_options)
        yaml_data = AOR_workflow.aor_workflow_setup(yaml_data, strata, automation_input)
    else:
        # Add new workflows here
        pass

    return yaml_data

def workflow_analysis(yaml_data, sm, analysis):
    """
    After the simulation has completed in openiam_cf.py, this function is called
    if a Workflow section is included. This function then checks for and runs
    a specific Workflow analysis.
    """
    if 'Type' in yaml_data['Workflow']:
        workflow_type = yaml_data['Workflow']['Type']

        if workflow_type == 'AoR':
            AOR_workflow.aor_workflow_analysis(yaml_data, sm, analysis)
        else:
            # Add future Workflow analysis types here.
            pass


