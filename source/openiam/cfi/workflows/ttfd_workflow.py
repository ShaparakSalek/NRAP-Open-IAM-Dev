
import logging

TTFD_AQUIFER_COMPONENT_OUTPUT = {
    'FutureGen2Aquifer': {'pH': ['pH_dx', 'pH_dy', 'pH_dz'],
                          'TDS': ['TDS_dx', 'TDS_dy', 'TDS_dz'],
                          'Dissolved_CO2': ['Dissolved_CO2_dx',
                                            'Dissolved_CO2_dy',
                                            'Dissolved_CO2_dz'],
                          'Pressure': ['Pressure_dx', 'Pressure_dy', 'Pressure_dz'],
                          'Temperature': ['Temperature_dx',
                                          'Temperature_dy',
                                          'Temperature_dz']},
    'FutureGen2AZMI': {'pH': ['pH_dx', 'pH_dy', 'pH_dz'],
                       'TDS': ['TDS_dx', 'TDS_dy', 'TDS_dz'],
                       'Dissolved_CO2': ['Dissolved_CO2_dx', 'Dissolved_CO2_dy',
                                         'Dissolved_CO2_dz'],
                       'Pressure': ['Pressure_dx', 'Pressure_dy', 'Pressure_dz'],
                       'Temperature': ['Temperature_dx',
                                       'Temperature_dy',
                                       'Temperature_dz']},
    'GenericAquifer': {'Dissolved_CO2': ['Dissolved_CO2_dr', 'Dissolved_CO2_dz'],
                       'Dissolved_salt': ['Dissolved_salt_dr', 'Dissolved_salt_dz']},
    'CarbonateAquifer': {'CarbonateAquifer': ['dx', 'dy']},
    'DeepAlluviumAquifer': {'pH': ['pH_dx', 'pH_dy', 'pH_dz'],
                            'TDS': ['TDS_dx', 'TDS_dy', 'TDS_dz'],
                            'Pressure': ['Pressure_dx', 'Pressure_dy', 'Pressure_dz']},
    'DeepAlluviumAquiferML': {'pH': ['pH_dx', 'pH_dy', 'pH_dz'],
                              'TDS': ['TDS_dx', 'TDS_dy', 'TDS_dz'],
                              'Pressure': ['Pressure_dx', 'Pressure_dy', 'Pressure_dz']}
    }

TTFD_AQUIFER_COMPONENT_PLUMES = {
    'FutureGen2Aquifer': ['pH', 'TDS', 'Dissolved_CO2', 'Pressure', 'Temperature'],
    'FutureGen2AZMI': ['pH', 'TDS', 'Dissolved_CO2', 'Pressure', 'Temperature'],
    'GenericAquifer': ['Dissolved_CO2', 'Dissolved_salt'],
    'CarbonateAquifer': ['CarbonateAquifer'],
    'DeepAlluviumAquifer': ['pH', 'TDS', 'Pressure'],
    'DeepAlluviumAquiferML': ['pH', 'TDS', 'Pressure']
    }

TTFD_DEFAULT_PLUME_TYPES = {
    'FutureGen2Aquifer': 'pH',
    'FutureGen2AZMI': 'pH',
    'GenericAquifer': 'Dissolved_CO2',
    'CarbonateAquifer': 'CarbonateAquifer',
    'DeepAlluviumAquifer': 'pH',
    'DeepAlluviumAquiferML': 'pH'
    }

TTFD_PLUME_TYPE_INCOMPATIBLE = ''.join([
    'The PlumeType entry provided for the TTFD workflow ({}) is incompatible ',
    'with the aquifer component type being used ({}). The simulation will used ',
    'the default plume type for that component type ({}).',
    ])

TTFD_PLUME_TYPE_MISSING = ''.join([
    'The TTFD workflow is being used, but the PlumeType entry was not present ',
    'in the .yaml file. For the aquifer component type being used ({}), the ',
    'default plume type of {} will be used.'])


def set_up_ttfd_workflow_plots(yaml_data, aquifer_component_type, figure_dpi=100,
                               plot_injection_sites=False, InjectionCoordx=None,
                               InjectionCoordy=None, fig_size=None, extension=None):
    """
    Sets up the plot entry dictionary for the TTFD workflow.
    """
    plumeType = yaml_data['Workflow']['Options'].get(
            'PlumeType', TTFD_DEFAULT_PLUME_TYPES[aquifer_component_type])

    compName = aquifer_component_type + '1'

    if not extension:
        extension= '.png'

    plotName = 'TTFD_Plot' + extension

    yaml_data['Plots'][plotName] = {'TTFD': {
        'PlumeType': plumeType, 'ComponentNameList': [compName],
        'PlotInjectionSites': plot_injection_sites,
        'FigureDPI': figure_dpi}}

    if InjectionCoordx is not None and InjectionCoordy is not None:
        yaml_data['Plots'][plotName]['TTFD'][
            'InjectionCoordx'] = InjectionCoordx
        yaml_data['Plots'][plotName]['TTFD'][
            'InjectionCoordy'] = InjectionCoordy

    if 'WriteDreamOutput' in yaml_data['Workflow']['Options']:
        dream_output = yaml_data['Workflow']['Options']['WriteDreamOutput']

        yaml_data['Plots'][plotName]['TTFD'][
            'WriteDreamOutput'] = dream_output

    if 'MonitoringLocations' in yaml_data['Workflow']['Options']:
        monitors = yaml_data['Workflow']['Options']['MonitoringLocations']

        yaml_data['Plots'][plotName]['TTFD'][
            'MonitoringLocations'] = monitors

    if 'NumZPointsWithinAquifers' in yaml_data['Workflow']['Options']:
        numZPointsAq = yaml_data['Workflow']['Options']['NumZPointsWithinAquifers']

        yaml_data['Plots'][plotName]['TTFD'][
            'NumZPointsWithinAquifers'] = numZPointsAq

    if 'NumZPointsWithinShales' in yaml_data['Workflow']['Options']:
        numZPointsSh = yaml_data['Workflow']['Options']['NumZPointsWithinShales']

        yaml_data['Plots'][plotName]['TTFD'][
            'NumZPointsWithinShales'] = numZPointsSh

    if 'SpecifyXandYLims' in yaml_data['Workflow']['Options']:
        lims = yaml_data['Workflow']['Options']['SpecifyXandYLims']

        yaml_data['Plots'][plotName]['TTFD'][
            'SpecifyXandYLims'] = lims

    if 'SpecifyXandYGridLims' in yaml_data['Workflow']['Options']:
        gridLims = yaml_data['Workflow']['Options']['SpecifyXandYGridLims']

        yaml_data['Plots'][plotName]['TTFD'][
            'SpecifyXandYGridLims'] = gridLims

    if 'xGridSpacing' in yaml_data['Workflow']['Options']:
        yaml_data['Plots'][plotName]['TTFD'][
            'xGridSpacing'] = yaml_data['Workflow']['Options']['xGridSpacing']

    if 'yGridSpacing' in yaml_data['Workflow']['Options']:
        yaml_data['Plots'][plotName]['TTFD'][
            'yGridSpacing'] = yaml_data['Workflow']['Options']['yGridSpacing']

    if fig_size:
        yaml_data['Plots'][plotName]['TTFD'][
            'FigureSize'] = fig_size

    if 'SaveCSVFiles' in yaml_data['Workflow']['Options']:
        yaml_data['Plots'][plotName]['TTFD'][
            'SaveCSVFiles'] = yaml_data['Workflow']['Options']['SaveCSVFiles']

    return yaml_data


def get_ttfd_aq_output_list(yaml_data, aquifer_component_type):
    """
    Returns the aquifer component outputs for the aquifer component type and
    plume type used for the TTFD workflow.
    """
    plume_type = yaml_data['Workflow']['Options'].get(
        'PlumeType', TTFD_DEFAULT_PLUME_TYPES[aquifer_component_type])

    output_list = TTFD_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][plume_type]

    return output_list


def set_up_ttfd_workflow(yaml_data, aquifer_component_type):
    """
    Checks the yaml dictionary for the options required for the TTFD workflow.
    """
    if 'PlumeType' in yaml_data['Workflow']['Options']:
        plumeType = yaml_data['Workflow']['Options']['PlumeType']

        if not plumeType in TTFD_AQUIFER_COMPONENT_PLUMES[aquifer_component_type]:
            logging.warning(TTFD_PLUME_TYPE_INCOMPATIBLE.format(
                aquifer_component_type, plumeType,
                TTFD_DEFAULT_PLUME_TYPES[aquifer_component_type]))

            plumeType = TTFD_DEFAULT_PLUME_TYPES[aquifer_component_type]

            yaml_data['Workflow']['Options']['PlumeType'] = plumeType

    else:
        plumeType = TTFD_DEFAULT_PLUME_TYPES[aquifer_component_type]

        yaml_data['Workflow']['Options']['PlumeType'] = plumeType

        logging.warning(TTFD_PLUME_TYPE_MISSING.format(aquifer_component_type, plumeType))

    return yaml_data
