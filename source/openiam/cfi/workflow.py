
import os
import sys
import logging

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)

try:
    from openiam import IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

import openiam.cfi.strata as iam_strata

# Certain default settings are stored in workflow_defaults.py
from openiam.cfi.workflow_defaults import (
    DEFAULT_RESERVOIR_COMP, DEFAULT_WELLBORE_COMP, DEFAULT_AQUIFER_COMP, 
    DEFAULT_FIGURE_DPI)
    
from openiam.cfi.workflow_utils import (
    get_parameter_names, get_automation_input, YAML_INPUT_WARNING_MSG)

from openiam.cfi.workflows.leakage_assessment_workflow import (
    set_up_leakage_assessment_workflow_plots)

from openiam.cfi.workflows.aor_workflow import (
    get_aor_loc_data, aor_crit_pressure_check, aor_workflow_analysis, 
    set_up_aor_workflow_plots, get_aor_aq_output_list, set_up_aor_workflow)

from openiam.cfi.workflows.ttfd_workflow import (
    set_up_ttfd_workflow_plots, get_ttfd_aq_output_list, set_up_ttfd_workflow)

# If new workflows are added, some sections of the code needs to be updated. I 
# used the term NEW_WORKFLOWS in those areas so they can be found more easily.
WORKFLOW_OPTIONS = ['LeakageAssessment', 'AoR', 'TTFD']

WORKFLOW_USE_RESERVOIR = {'LeakageAssessment': True, 'AoR': True, 
                          'TTFD': True}

WORKFLOW_USE_WELLBORE = {'LeakageAssessment': True, 'AoR': True, 
                         'TTFD': True}

WORKFLOW_USE_AQUIFER = {'LeakageAssessment': False, 'AoR': True, 
                        'TTFD': True}

# Specifies whether specific options need to be handled for the workflow
ADDITIONAL_SETUP_REQUIRED = {'LeakageAssessment': False, 'AoR': True, 
                             'TTFD': True}

# Specifies whether a workflow requires location input
LOC_INPUT_REQUIRED = {'LeakageAssessment': True, 'AoR': False, 
                      'TTFD': True}

# The reservoir and wellbore outputs for each component type.
WORKFLOW_RESERVOIR_OUTPUT = {
    'LeakageAssessment': ['pressure', 'CO2saturation'], 
    'AoR': ['pressure', 'CO2saturation'], 
    'TTFD': ['pressure', 'CO2saturation']
    }

WORKFLOW_WELLBORE_OUTPUT = {
    'LeakageAssessment': ['brine_aquifer{}', 'CO2_aquifer{}'], 
    'AoR': ['brine_aquifer{}', 'CO2_aquifer{}'], 
    'TTFD': ['brine_aquifer{}', 'CO2_aquifer{}']
    }

LOC_INPUT_ERR_MSG = ''.join([
    'No Location data was provided under the WellboreOptions entry ', 
    '(contained in the Workflow: Options section of the .yaml file). ', 
    'The Workflow type was {}, however, and the {} Workflow ',
    'requires the user to provide Location data for the wellbores. ', 
    'This simulation cannot proceed.'])

DEFAULT_LUTR_FILE_DIRECTORY = os.path.join(IAM_DIR, 'source', 'components',
                                           'reservoir', 'lookuptables',
                                           'FutureGen2', '1008_sims')

DEFAULT_LUTR_ENTRIES = {
    'FileDirectory': DEFAULT_LUTR_FILE_DIRECTORY,
    'TimeFile': 'time_points.csv',
    'ParameterFilename': 'parameters_and_filenames_trunc.csv',
    'Interpolation2D': True, 'Parameters': {'index': 1}}


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
        if yaml_data['Workflow']['Type'] not in WORKFLOW_OPTIONS:
            warning_msg = ''.join([
                'A Workflow section was included in the .yaml file, but the Workflow ',
                'Type provided (', yaml_data['Workflow']['Type'], ') was not ',
                'recognized as one of the available options (', str(WORKFLOW_OPTIONS),
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

    yaml_data = workflow_setup(yaml_data, strata)

    return yaml_data


def workflow_setup(yaml_data, strata):
    """
    Sets up the components and plots required for the Workflow.
    """
    types_strata_pars = iam_strata.get_comp_types_strata_pars()
    types_strata_obs = iam_strata.get_comp_types_strata_obs()
    
    workflow_type = yaml_data['Workflow']['Type']
    
    automation_input = get_automation_input(yaml_data)

    res_component_type = yaml_data['Workflow']['Options'].get(
        'ReservoirComponentType', DEFAULT_RESERVOIR_COMP)

    well_component_type = yaml_data['Workflow']['Options'].get(
        'WellboreComponentType', DEFAULT_WELLBORE_COMP)

    aquifer_component_type = yaml_data['Workflow']['Options'].get(
        'AquiferComponentType', DEFAULT_AQUIFER_COMP)
    
    strata_type = iam_strata.get_strata_type_from_yaml(yaml_data)
    
    if strata_type in types_strata_pars:
        numberOfShaleLayers = strata.get_num_shale_layers()
    elif strata_type in types_strata_obs:
        numberOfShaleLayers = 3
        
        strata_data = yaml_data[strata_type]
        
        if 'Parameters' in strata_data:
            if 'numberOfShaleLayers' in strata_data['Parameters']:
                numberOfShaleLayersInput = strata_data['Parameters']['numberOfShaleLayers']
                if isinstance(numberOfShaleLayersInput, dict):
                    if 'value' in numberOfShaleLayersInput:
                        numberOfShaleLayers = numberOfShaleLayersInput['value']
                else:
                    numberOfShaleLayers = numberOfShaleLayersInput

    default_aquifer_name = 'aquifer{}'.format(str(numberOfShaleLayers - 1))

    aquifer_name = yaml_data['Workflow']['Options'].get(
        'AquiferName', default_aquifer_name)
    
    # Check if AquiferName was included under AquiferOptions - if so, use that 
    # input and log a warning message.
    if 'AquiferOptions' in yaml_data['Workflow']['Options']:
        if 'AquiferName' in yaml_data['Workflow']['Options']['AquiferOptions']:
            aq_name_2 = yaml_data['Workflow']['Options']['AquiferOptions']
            
            if aquifer_name != aq_name_2:
                warning_msg_pt1 = ''.join([
                    'The AquiferName input was provided under AquiferOptions ',
                    'entry (AquiferName: {}). '.format(aq_name_2)])
                
                if 'AquiferName' in yaml_data['Workflow']['Options']:
                    warning_msg_pt2 = ''.join([
                        'The AquiferName input was also provided under the Options ',
                        'entry, however (AquiferName: {}). '.format(aquifer_name)])
                else:
                    warning_msg_pt2 = ''.join([
                        'The AquiferName input was not provided under the Options ',
                        'entry. '])
                
                warning_msg_pt3 = ''.join([
                    'The AquiferName input provided under AquiferOptions ',
                    'will be used (AquiferName: {}). '.format(aq_name_2)])
                
                warning_msg = warning_msg_pt1 + warning_msg_pt2 + warning_msg_pt3
                
                logging.warning(warning_msg)
                
                aquifer_name = aq_name_2
    
    # This section is meant for handling special requirements for the workflows 
    # (e.g., required options). Add cases for NEW_WORKFLOWS here, if neccessary.
    if ADDITIONAL_SETUP_REQUIRED[workflow_type]:
        if workflow_type == 'AoR':
            yaml_data = set_up_aor_workflow(yaml_data)
            
        elif workflow_type == 'TTFD':
            yaml_data = set_up_ttfd_workflow(yaml_data, aquifer_component_type)
    
    # Now, handle the components
    if 'Components' not in yaml_data['ModelParams']:
        yaml_data['ModelParams']['Components'] = []

    elif not isinstance(yaml_data['ModelParams']['Components'], list):
        comps = yaml_data['ModelParams']['Components']

        yaml_data['ModelParams']['Components'] = []
        for comp in comps:
            yaml_data['ModelParams']['Components'].append(comp)
    
    # Add the component entries for the workflow
    if automation_input['ResComp'] and WORKFLOW_USE_RESERVOIR[workflow_type]:
        yaml_data, res_component_name = add_res_component_entries(
            yaml_data, res_component_type, workflow_type)

        yaml_data['ModelParams']['Components'].append(res_component_name)

    if automation_input['WellComp'] and WORKFLOW_USE_WELLBORE[workflow_type]:
        yaml_data, well_component_name = add_well_component_entries(
            yaml_data, well_component_type, aquifer_name, res_component_name, 
            workflow_type)

        yaml_data['ModelParams']['Components'].append(well_component_name)

    if automation_input['AqComp'] and WORKFLOW_USE_AQUIFER[workflow_type]:
        yaml_data, aq_component_name = add_aq_component_entries(
            yaml_data, aquifer_component_type, aquifer_name, well_component_name, 
            workflow_type)

        yaml_data['ModelParams']['Components'].append(aq_component_name)
    
    # Add the plot entries for the workflow
    if automation_input['Plots']:
        yaml_data = set_up_plots_section(yaml_data, well_component_type, 
                                         aquifer_component_type, 
                                         workflow_type, aquifer_name)

    return yaml_data


def workflow_analysis(yaml_data, sm, analysis):
    """
    After the simulation has completed in openiam_cf.py, this function is called
    if a Workflow section is included. This function then checks for and runs
    a specific Workflow analysis. Some, like TTFD, do not need to be run separately, 
    as the TTFD plot is already run after the simulation (i.e., it is in the 
    Plots section of the .yaml file). This function is for analyses that need 
    to be run after everything else has finished.
    """
    if 'Type' in yaml_data['Workflow']:
        # Add cases for NEW_WORKFLOWS here, if necessary.
        if yaml_data['Workflow']['Type'] == 'AoR':
            aor_workflow_analysis(yaml_data, sm, analysis)


def add_res_component_entries(yaml_data, res_component_type, workflow_type):
    """
    This automatically sets up the reservoir component entry for the Workflow,
    if AutomateResCompSetup is set to True.
    """
    compName = res_component_type + '1'

    output_list = WORKFLOW_RESERVOIR_OUTPUT[workflow_type] 

    yaml_data[compName] = {'Type': res_component_type, 'Outputs': output_list}

    if res_component_type == 'LookupTableReservoir':
        lutr_entries = DEFAULT_LUTR_ENTRIES

        if 'ReservoirOptions' in yaml_data['Workflow']['Options']:
            lutr_entries.update(
                yaml_data['Workflow']['Options']['ReservoirOptions'])

        for entry in list(lutr_entries.keys()):
            yaml_data[compName][entry] = lutr_entries[entry]

    if 'ReservoirOptions' in yaml_data['Workflow']['Options'] and \
            res_component_type != 'LookupTableReservoir':
        if 'InjectionWell' in yaml_data['Workflow']['Options']['ReservoirOptions']:
            injection_x = yaml_data['Workflow']['Options']['ReservoirOptions'][
                'InjectionWell']['coordx']
            injection_y = yaml_data['Workflow']['Options']['ReservoirOptions'][
                'InjectionWell']['coordy']

            yaml_data[compName]['InjectionWell'] = {'coordx': injection_x,
                                                    'coordy': injection_y}

        if 'Parameters' in yaml_data['Workflow']['Options']['ReservoirOptions']:
            par_names = get_parameter_names(res_component_type)

            first_time_check = True
            for par in par_names:
                if par in yaml_data['Workflow']['Options']['ReservoirOptions'][
                        'Parameters']:
                    if first_time_check:
                        yaml_data[compName]['Parameters'] = {}
                        first_time_check = False

                    yaml_data[compName]['Parameters'][par] = yaml_data['Workflow'][
                        'Options']['ReservoirOptions']['Parameters'][par]

    return yaml_data, compName


def add_well_component_entries(yaml_data, well_component_type, aquifer_name,
                               res_component_name, workflow_type):
    """
    This automatically sets up the component entries for the workflow, if
    AutomateWellCompSetup is set to True.
    """
    compName = well_component_type + '1'

    aquiferNum = aquifer_name[7:]
    
    workflow_output = WORKFLOW_WELLBORE_OUTPUT[workflow_type]

    output_list = [val.format(aquiferNum) for val in workflow_output]
    
    controls = {'critPressureApproach': True, 'enforceCritPressure': False}
    
    loc_data_check = False
    if 'WellboreOptions' in yaml_data['Workflow']['Options']:
        if 'Controls' in yaml_data['Workflow']['Options']['WellboreOptions']:
            controls.update(yaml_data['Workflow']['Options']['WellboreOptions']['Controls'])

        if 'Locations' in yaml_data['Workflow']['Options']['WellboreOptions']:
            loc_data = yaml_data['Workflow']['Options']['WellboreOptions']['Locations']
            loc_data_check = True
    
    # Add more cases for NEW_WORKFLOWS here. Some workflows require Locations 
    # data (e.g., LeakageAssessment), while some do not (AoR).
    if not loc_data_check:
        if LOC_INPUT_REQUIRED[workflow_type]:
            # Locations input is required but none was given, resulting in an error
            logging.error(LOC_INPUT_ERR_MSG.format(workflow_type, workflow_type))
        
            raise KeyError(LOC_INPUT_ERR_MSG.format(workflow_type, workflow_type))
        
        elif not LOC_INPUT_REQUIRED[workflow_type]:
            # Locations input was not given, but it is not required for the workflow
            if workflow_type == 'AoR':
                loc_data = get_aor_loc_data(yaml_data, res_component_name, loc_data_check)
    
    if workflow_type == 'AoR':
        # This function checks if CriticalPressureMPa was set to a specific 
        # pressure. If it was AND an OpenWellbore component is being used, the 
        # wellbore is made to use that critical pressure as the critPressure 
        # parameter. The controls dictionary also has to be updated.
        yaml_data, controls = aor_crit_pressure_check(
            yaml_data, well_component_type, controls)
    
    yaml_data[compName] = {'Type': well_component_type, 'LeakTo': aquifer_name,
                           'Connection': res_component_name, 'Locations': loc_data,
                           'Outputs': output_list}
    
    if well_component_type == 'OpenWellbore':
        yaml_data[compName]['Controls'] = controls
    
    if 'WellboreOptions' in yaml_data['Workflow']['Options']:
        if 'Controls' in yaml_data['Workflow']['Options']['WellboreOptions'] \
            and well_component_type != 'OpenWellbore':
                yaml_data[compName]['Controls'] = controls
        
        if 'LeakTo' in yaml_data['Workflow']['Options']['WellboreOptions']:
            yaml_data[compName]['LeakTo'] = yaml_data['Workflow']['Options'][
                'WellboreOptions']['LeakTo']
        
        if 'ThiefZone' in yaml_data['Workflow']['Options']['WellboreOptions']:
            yaml_data[compName]['ThiefZone'] = yaml_data['Workflow']['Options'][
                'WellboreOptions']['ThiefZone']
        
        if 'Parameters' in yaml_data['Workflow']['Options']['WellboreOptions']:
            par_names = get_parameter_names(well_component_type)

            first_time_check = True
            for par in par_names:
                if par in yaml_data['Workflow']['Options']['WellboreOptions'][
                        'Parameters']:
                    if first_time_check:
                        yaml_data[compName]['Parameters'] = {}
                        first_time_check = False

                    yaml_data[compName]['Parameters'][par] = yaml_data['Workflow'][
                        'Options']['WellboreOptions']['Parameters'][par]

            # critPressure will not be in the default_pars.keys(), which is used
            # in get_parameter_names(). So check that parameter separately.
            if 'critPressure' in yaml_data['Workflow']['Options']['WellboreOptions'][
                    'Parameters']:
                yaml_data[compName]['Parameters']['critPressure'] = yaml_data[
                    'Workflow']['Options']['WellboreOptions']['Parameters'][
                        'critPressure']

    return yaml_data, compName


def add_aq_component_entries(yaml_data, aquifer_component_type, aquifer_name,
                             well_component_name, workflow_type):
    """
    This automatically sets up the component entries for the workflow, if
    AutomateAqCompSetup is set to True.
    """
    compName = aquifer_component_type + '1'
    
    # Add cases for NEW_WORKFLOWS here
    if workflow_type == 'TTFD':
        output_list = get_ttfd_aq_output_list(yaml_data, aquifer_component_type)
        
    elif workflow_type == 'AoR':
        output_list = get_aor_aq_output_list(yaml_data, aquifer_component_type)
        
    yaml_data[compName] = {'Type': aquifer_component_type, 'AquiferName': aquifer_name,
                           'Connection': well_component_name, 'Outputs': output_list}

    if 'AquiferOptions' in yaml_data['Workflow']['Options']:
        if 'Parameters' in yaml_data['Workflow']['Options']['AquiferOptions']:
            par_names = get_parameter_names(aquifer_component_type)

            first_time_check = True
            for par in par_names:
                if par in yaml_data['Workflow']['Options']['AquiferOptions'][
                        'Parameters']:
                    if first_time_check:
                        yaml_data[compName]['Parameters'] = {}
                        first_time_check = False

                    yaml_data[compName]['Parameters'][par] = yaml_data['Workflow'][
                        'Options']['AquiferOptions']['Parameters'][par]

    return yaml_data, compName


def set_up_plots_section(yaml_data, well_component_type, aquifer_component_type, 
                         workflow_type, aquifer_name):
    """
    This automatically sets up the Plots section of the .yaml file for the 
    Workflow, if AutomatePlotsSetup is set to True.
    """
    # These options are used by multiple workflows, so they are obtained here 
    # to avoid any repetition within separate workflow files.
    figure_dpi = yaml_data['Workflow']['Options'].get('FigureDPI', DEFAULT_FIGURE_DPI)
    
    extension = ''
    if 'FigureName' in yaml_data['Workflow']['Options']:
        fig_name = yaml_data['Workflow']['Options']['FigureName']
        if '.' in fig_name:
            extension = fig_name[fig_name.index('.'):]
    
    fig_size = None
    if 'FigureSize' in yaml_data['Workflow']['Options']:
        fig_size = yaml_data['Workflow']['Options']['FigureSize']
    
    plot_injection_sites = False
    if 'PlotInjectionSites' in yaml_data['Workflow']['Options']:
        if isinstance(yaml_data['Workflow']['Options']['PlotInjectionSites'], bool):
            plot_injection_sites = yaml_data['Workflow']['Options']['PlotInjectionSites']
        else:
            warning_msg = YAML_INPUT_WARNING_MSG.format(
                'PlotInjectionSites', 'boolean', 'False')
            logging.warning(warning_msg)
    
    InjectionCoordx = None
    InjectionCoordy = None
    if 'InjectionCoordx' in yaml_data['Workflow']['Options'] and \
            'InjectionCoordy' in yaml_data['Workflow']['Options']:
        InjectionCoordx = yaml_data['Workflow']['Options']['InjectionCoordx']
        InjectionCoordy = yaml_data['Workflow']['Options']['InjectionCoordy']

    if 'Plots' not in yaml_data:
        yaml_data['Plots'] = {}
    
    # Add cases for NEW_WORKFLOWS here
    if workflow_type == 'LeakageAssessment':
        yaml_data = set_up_leakage_assessment_workflow_plots(
            yaml_data, aquifer_name, figure_dpi=figure_dpi, fig_size=fig_size, 
            extension=extension)
        
    elif workflow_type == 'AoR':
        yaml_data = set_up_aor_workflow_plots(
            yaml_data, well_component_type, aquifer_component_type, 
            figure_dpi=figure_dpi, plot_injection_sites=plot_injection_sites, 
            InjectionCoordx=InjectionCoordx, InjectionCoordy=InjectionCoordy, 
            fig_size=fig_size, extension=extension)
        
    elif workflow_type == 'TTFD':
        yaml_data = set_up_ttfd_workflow_plots(
            yaml_data, aquifer_component_type, figure_dpi=figure_dpi, 
            plot_injection_sites=plot_injection_sites, 
            InjectionCoordx=InjectionCoordx, InjectionCoordy=InjectionCoordy, 
            fig_size=fig_size, extension=extension)
    
    return yaml_data
