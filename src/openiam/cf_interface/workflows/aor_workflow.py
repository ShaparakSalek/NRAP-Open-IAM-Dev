
import os
import sys
import logging
import csv
import numpy as np
import pandas as pd

try:
    from openiam.components.iam_base_classes import IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.components.open_wellbore_component import OpenWellbore

from openiam.cf_interface.strata import (get_strata_type_from_yaml,
                                         get_comp_types_strata_pars,
                                         get_comp_types_strata_obs)

import openiam.visualization.area_of_review as AoR
from openiam.visualization.area_of_review import CSV_FILE_NAME_TAGS as AOR_CSV_FILE_NAME_TAGS
from openiam.visualization.area_of_review import CSV_FILE_COLUMNS as AOR_CSV_FILE_COLUMNS


from openiam.cf_interface.workflow_defaults import DEFAULT_AQUIFER_COMP

DEFAULT_CRIT_PRESSURE_SETTING = 'Calculated'

AOR_PLOT_NAMES = [
    {'Pressure_Plot': 'pressure'},
    {'CO2_Saturation_Plot': 'CO2saturation'},
    {'Aq_CO2_Impact_Plot': '{AqCO2Metric}'},
    {'Aq_Brine_Impact_Plot': '{AqSaltMetric}'}]

AOR_AQUIFER_COMPONENTS = ['FutureGen2Aquifer', 'FutureGen2AZMI',
                          'GenericAquifer', 'CarbonateAquifer',
                          'DeepAlluviumAquifer', 'DeepAlluviumAquiferML']

AOR_CRIT_PRESSURE_COLUMN = 'Critical Pressure [MPa]'

AOR_AQUIFER_COMPONENT_OUTPUT = {
    'FutureGen2Aquifer': ['pH_volume', 'TDS_volume'],
    'FutureGen2AZMI': ['pH_volume', 'TDS_volume'],
    'GenericAquifer': ['Dissolved_CO2_volume', 'Dissolved_salt_volume'],
    'CarbonateAquifer': ['pH_volume', 'TDS_volume'],
    'DeepAlluviumAquifer': ['pH_volume', 'TDS_volume'],
    'DeepAlluviumAquiferML': ['pH_volume', 'TDS_volume']
    }

# Wellbore component types that have a brine density parameter (e.g., brineDensity).
# If this component type is used AND a BrineDensity input is given under
# yaml_data['Workflow']['Options'], the the BrineDensity input will be removed.
WELL_COMPS_WITH_BRINE_DENSITY = ['OpenWellbore', 'MultisegmentedWellbore']

# These are the default minimum spacings to use when thinning the point densities
# from LookupTableReservoir files.
DEFAULT_MIN_X_SPACING = 20000
DEFAULT_MIN_Y_SPACING = 20000

# Default grid limits when using the AoR workflow and a reservoir component
# that is not a LookupTableReservoir.
DEFAULT_AOR_XMIN = -50000
DEFAULT_AOR_XMAX = 50000
DEFAULT_AOR_XSIZE = 6

DEFAULT_AOR_YMIN = -50000
DEFAULT_AOR_YMAX = 50000
DEFAULT_AOR_YSIZE = 6


def get_aor_loc_data(yaml_data, res_component_name, loc_data_check):
    """
    When the AoR workflow is used and the user does not provide location input,
    this function returns the loc_data dictionary for add_well_component_entries()
    in workflow.py. The location data are either obtained from the LookupTableReservoir
    .csv files or taken as a grid.
    """
    if yaml_data[res_component_name]['Type'] == 'LookupTableReservoir' and not loc_data_check:
        file_dir = yaml_data[res_component_name]['FileDirectory']

        par_file_name = yaml_data[res_component_name]['ParameterFilename']

        index = int(yaml_data[res_component_name]['Parameters']['index'])

        data = pd.read_csv(os.path.join(IAM_DIR, file_dir, par_file_name))

        filename = data['filename'][index]

        file_name_dir = os.path.join(IAM_DIR, file_dir, filename)

        read_z_values = False

        data = pd.read_csv(file_name_dir)

        if 'z' in data:
            read_z_values = True

        del data

        # If no location data are given for the OpenWelbores, use the LUTR file itself by default
        thin_point_density = yaml_data['Workflow']['Options'].get(
            'thin_point_density', True)
        min_x_spacing = yaml_data['Workflow']['Options'].get(
            'min_x_spacing', DEFAULT_MIN_X_SPACING)
        min_y_spacing = yaml_data['Workflow']['Options'].get(
            'min_y_spacing', DEFAULT_MIN_Y_SPACING)

        loc_data = {'file': file_name_dir, 'read_z_values': read_z_values,
                    'thin_point_density': thin_point_density,
                    'min_x_spacing': min_x_spacing, 'min_y_spacing': min_y_spacing}
        loc_data_check = True

    # If no wellbore locations were given and the simulation does not use a
    # LookupTableReservoir (from which to draw wellbore locations), use default values
    if not loc_data_check:
        xmin = DEFAULT_AOR_XMIN
        xmax = DEFAULT_AOR_XMAX
        xsize = DEFAULT_AOR_XSIZE
        ymin = DEFAULT_AOR_YMIN
        ymax = DEFAULT_AOR_YMAX
        ysize = DEFAULT_AOR_YSIZE
        loc_data = {'grid': {'xmin': xmin, 'xmax': xmax, 'xsize': xsize,
                             'ymin': ymin, 'ymax': ymax, 'ysize': ysize}}

    return loc_data


def get_points_in_aor(yaml_data, sm=None, time_index=None):
    """
    Checks the AoR plot type results for multiple metrics and determines an
    AoR that reflects all metrics considered. Pressure results are only considered
    if a critical pressure is provided.
    """
    # These lists indicate the stratigraphy component types that offer thicknesses
    # and depths as parameters or as observations.
    types_strata_pars = get_comp_types_strata_pars()
    types_strata_obs = get_comp_types_strata_obs()

    output_dir = yaml_data['ModelParams']['OutputDirectory']

    aquifer_component_type = yaml_data['Workflow']['Options'].get(
        'AquiferComponentType', DEFAULT_AQUIFER_COMP)

    # Get the stratigraphy information from the .yaml file
    strata_type = get_strata_type_from_yaml(yaml_data)

    AoR_included_x_km = []
    AoR_included_y_km = []

    if time_index is None:
        AoR_filename = 'AoR_{}.csv'
    else:
        AoR_filename = 'AoR_{}' + '_tIndex_{}.csv'.format(time_index)

    pressure_included = True
    # Reservoir pressures, only use if a critical pressure was provided
    if 'CriticalPressureMPa' in yaml_data['Workflow']['Options']:
        # Reservoir pressures
        pressure_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS['pressure'])

        AoR_pressure_results = pd.read_csv(os.path.join(IAM_DIR, output_dir,
                                                        'csv_files', pressure_file))

        pressure_x_km = AoR_pressure_results['x (km)'].values
        pressure_y_km = AoR_pressure_results['y (km)'].values

        pressure_vals_MPa = AoR_pressure_results[AOR_CSV_FILE_COLUMNS['pressure']].values

        try:
            if strata_type in types_strata_pars:
                critPressureVal = np.min(AoR_pressure_results[AOR_CRIT_PRESSURE_COLUMN])
            elif strata_type in types_strata_obs:
                critPressureVal = AoR_pressure_results[AOR_CRIT_PRESSURE_COLUMN]
        except:
            if yaml_data['Workflow']['Options']['CriticalPressureMPa'] == 'Calculated':
                # returned in MPa
                critPressureVal = get_crit_pressure_aor_analysis(
                    len(pressure_x_km), yaml_data, sm)
            else:
                critPressureVal = float(yaml_data['Workflow']['Options']['CriticalPressureMPa'])

        if strata_type in types_strata_pars:
            pressure_x_km_AoR = pressure_x_km[pressure_vals_MPa >= critPressureVal]
            pressure_y_km_AoR = pressure_y_km[pressure_vals_MPa >= critPressureVal]

            if len(pressure_x_km_AoR) > 0:
                AoR_included_x_km += pressure_x_km_AoR.tolist()
                AoR_included_y_km += pressure_y_km_AoR.tolist()

        elif strata_type in types_strata_obs:
            for loc_ref in range(len(pressure_x_km)):
                if pressure_vals_MPa[loc_ref] >= critPressureVal[loc_ref]:
                    AoR_included_x_km.append(pressure_x_km[loc_ref])
                    AoR_included_y_km.append(pressure_y_km[loc_ref])

    else:
        pressure_x_km = None
        pressure_y_km = None

        pressure_included = False

    # Reservoir CO2 saturations
    CO2saturation_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS['CO2saturation'])

    AoR_CO2saturation_results = pd.read_csv(os.path.join(
        IAM_DIR, output_dir, 'csv_files', CO2saturation_file))

    CO2saturation_x_km = AoR_CO2saturation_results['x (km)'].values
    CO2saturation_y_km = AoR_CO2saturation_results['y (km)'].values

    CO2saturation_vals = AoR_CO2saturation_results[AOR_CSV_FILE_COLUMNS['CO2saturation']].values

    CO2saturation_x_km_AoR = CO2saturation_x_km[CO2saturation_vals > 0]
    CO2saturation_y_km_AoR = CO2saturation_y_km[CO2saturation_vals > 0]

    if len(CO2saturation_x_km_AoR) > 0:
        AoR_included_x_km += CO2saturation_x_km_AoR.tolist()
        AoR_included_y_km += CO2saturation_y_km_AoR.tolist()

    # pH or dissolved CO2 contaminant plume volumes
    CO2impact_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][0]])

    AoR_CO2impact_results = pd.read_csv(os.path.join(IAM_DIR, output_dir,
                                                     'csv_files', CO2impact_file))

    CO2impact_x_km = AoR_CO2impact_results['x (km)'].values
    CO2impact_y_km = AoR_CO2impact_results['y (km)'].values

    CO2impact_vals_m3 = AoR_CO2impact_results[AOR_CSV_FILE_COLUMNS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][0]]].values

    CO2impact_x_km_AoR = CO2impact_x_km[CO2impact_vals_m3 > 0]
    CO2impact_y_km_AoR = CO2impact_y_km[CO2impact_vals_m3 > 0]

    if len(CO2impact_x_km_AoR) > 0:
        AoR_included_x_km += CO2impact_x_km_AoR.tolist()
        AoR_included_y_km += CO2impact_y_km_AoR.tolist()

    # TDS or dissolved salt contaminant plume volumes
    brineimpact_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][1]])

    AoR_brineimpact_results = pd.read_csv(os.path.join(IAM_DIR, output_dir,
                                                       'csv_files', brineimpact_file))

    brineimpact_x_km = AoR_brineimpact_results['x (km)'].values
    brineimpact_y_km = AoR_brineimpact_results['y (km)'].values

    brineimpact_vals_m3 = AoR_brineimpact_results[AOR_CSV_FILE_COLUMNS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][1]]].values

    brineimpact_x_km_AoR = brineimpact_x_km[brineimpact_vals_m3 > 0]
    brineimpact_y_km_AoR = brineimpact_y_km[brineimpact_vals_m3 > 0]

    if len(brineimpact_x_km_AoR) > 0:
        AoR_included_x_km += brineimpact_x_km_AoR.tolist()
        AoR_included_y_km += brineimpact_y_km_AoR.tolist()

    # Remove redundant entries
    AoR_included_x_km, AoR_included_y_km = remove_redundant_points(
        AoR_included_x_km, AoR_included_y_km)

    # Assemble the points that are not in the AoR.
    # To do that, evaluate all points considered. The x and y points for all
    # result types should be the same, this function just verifies that.
    All_x_points_km, All_y_points_km = verify_aor_points(
        pressure_x_km, pressure_y_km, CO2saturation_x_km, CO2saturation_y_km,
        CO2impact_x_km, CO2impact_y_km, brineimpact_x_km, brineimpact_y_km)

    AoR_point_included = np.zeros(len(All_x_points_km))

    for loc_ref, (x_loc, y_loc) in enumerate(zip(All_x_points_km, All_y_points_km)):
        AoR_point_temp = (x_loc, y_loc)

        point_included = False
        for (x_loc_in, y_loc_in) in zip(AoR_included_x_km, AoR_included_y_km):
            AoR_included_loc_temp = (x_loc_in, y_loc_in)

            if AoR_point_temp == AoR_included_loc_temp:
                point_included = True

        if point_included:
            AoR_point_included[loc_ref] = 1

    results_formatted = np.empty(((len(All_x_points_km) + 1), 3), dtype=list)

    results_formatted[0, 0] = 'Considered AoR Point, x (km)'
    results_formatted[0, 1] = 'Considered AoR Point, y (km)'
    results_formatted[0, 2] = 'Point Included in AoR'

    results_formatted[1:, 0] = All_x_points_km
    results_formatted[1:, 1] = All_y_points_km
    results_formatted[1:, 2] = AoR_point_included

    if time_index is not None:
        filename = os.path.join(output_dir, 'csv_files',
                                'AoR_Workflow_Output_tIndex{}.csv'.format(time_index))
    else:
        filename = os.path.join(output_dir, 'csv_files', 'AoR_Workflow_Output.csv')

    # Save the workflow outputs
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for row_ref in range(len(All_x_points_km) + 1):
            writer.writerow(results_formatted[row_ref, :])

    f.close()

    return All_x_points_km, All_y_points_km, AoR_point_included, pressure_included


def remove_redundant_points(AoR_x_km, AoR_y_km):
    """
    This function takes lists of x and y points and returns the lists after
    removing redundant combinations of x and y.
    """
    AoR_x_km_edit = []
    AoR_y_km_edit = []

    AoR_locs = []
    for (xLoc, yLoc) in zip(AoR_x_km, AoR_y_km):
        AoR_loc_temp = (xLoc, yLoc)

        if AoR_loc_temp not in AoR_locs:
            AoR_x_km_edit.append(xLoc)
            AoR_y_km_edit.append(yLoc)

            AoR_locs.append(AoR_loc_temp)

    AoR_x_km = AoR_x_km_edit
    AoR_y_km = AoR_y_km_edit

    return AoR_x_km, AoR_y_km


def verify_aor_points(pressure_x_km, pressure_y_km, CO2saturation_x_km,
                      CO2saturation_y_km, CO2impact_x_km, CO2impact_y_km,
                      brineimpact_x_km, brineimpact_y_km):
    """
    The x and y locations used in each of the plot types should be the same:
    this function ensures that they are indeed the same. If not, it produces an error.
    """
    x_all_same = True
    y_all_same = True

    metrics_string = '{}CO2 saturations{} and aquifer plume volumes from CO2 and brine Leakage'
    if pressure_x_km is None and pressure_y_km is None:
        x_all_same = (CO2saturation_x_km.all() == CO2impact_x_km.all() and
                      CO2impact_x_km.all() == brineimpact_x_km.all())

        y_all_same = (CO2saturation_y_km.all() == CO2impact_y_km.all() and
                      CO2impact_y_km.all() == brineimpact_y_km.all())

        metrics_string = metrics_string.format('', '')
    else:
        x_all_same = (pressure_x_km.all() == CO2saturation_x_km.all() and
                      CO2saturation_x_km.all() == CO2impact_x_km.all() and
                      CO2impact_x_km.all() == brineimpact_x_km.all())

        y_all_same = (pressure_y_km.all() == CO2saturation_y_km.all() and
                      CO2saturation_y_km.all() == CO2impact_y_km.all() and
                      CO2impact_y_km.all() == brineimpact_y_km.all())

        metrics_string = metrics_string.format('pressures, ', ',')

    if x_all_same and y_all_same:
        All_x_points_km = brineimpact_x_km
        All_y_points_km = brineimpact_y_km
    else:
        error_msg = ''.join([
            'The x and y points from the AoR .csv files for ', metrics_string,
            'were not identical. The x and y points for each AoR plot type ',
            'should be the same, if they are from one simulation with fixed ',
            'wellbore locations. Check your input for the simulation. The AoR ',
            'Workflow cannot be finished.'])
        logging.error(error_msg)

    return All_x_points_km, All_y_points_km


def get_crit_pressure_aor_analysis(num_pressure_points, yaml_data, sm):
    """
    Function that calculates the critical pressures (in MPa) for OpenWellbore
    components. When using spatially variable stratigraphy, the critical pressures
    are returned in an array where the value in each row corresponds with the
    pressure_x_km and pressure_y_km values in the same row.
    """
    types_strata_pars = get_comp_types_strata_pars()
    types_strata_obs = get_comp_types_strata_obs()

    # Get the stratigrapy information from the .yaml file
    strata_type = get_strata_type_from_yaml(yaml_data)

    if strata_type in types_strata_pars:
        critPressure = None
    elif strata_type in types_strata_obs:
        critPressure = np.zeros((num_pressure_points, 1))

    components = list(sm.component_models.values())

    for output_component in components:

        if isinstance(output_component, OpenWellbore):
            if strata_type in types_strata_pars and critPressure is None:
                # If using uniform stratigraphy, only do this once
                critPressureVal = AoR.get_crit_pressure(
                    output_component, sm=sm, yaml_data=yaml_data)

                critPressure = critPressureVal / 1.0e+6

            elif strata_type in types_strata_obs:
                critPressureVal = AoR.get_crit_pressure(
                    output_component, sm=sm, yaml_data=yaml_data)

                # Find the location index
                loc_ref = int(output_component.name.split('_')[-1])

                critPressure[loc_ref] = critPressureVal / 1.0e+6

    return critPressure


def aor_workflow_analysis(yaml_data, sm, analysis):
    """
    Evaluates the analysis of the data saved by each AoR plot entry and delineates
    an overall AoR that relects all metrics considered: pressure, CO2saturation,
    and both CO2 and brine impacts on the aquifer considered (i.e., plume volumes).
    """

    TimeList_option = False
    if 'TimeList' in yaml_data['Workflow']['Options']:
        if isinstance(yaml_data['Workflow']['Options']['TimeList'], list) or \
                yaml_data['Workflow']['Options']['TimeList'] == 'All':
            TimeList_option = True

            TimeList = yaml_data['Workflow']['Options']['TimeList']

            time = sm.time_array / 365.25
            if TimeList == 'All':
                time_index_list = range(len(time))
            else:
                time_index_list = AoR.get_t_indices(TimeList, time)

    if TimeList_option:
        for tIndex in time_index_list:
            All_x_points_km, All_y_points_km, AoR_point_included, pressure_included = \
                get_points_in_aor(yaml_data, sm=sm, time_index=tIndex)

            AoR.plot_aor_workflow_results(
                yaml_data, sm, All_x_points_km, All_y_points_km,
                AoR_point_included, time_index=tIndex, analysis=analysis)
    else:
        All_x_points_km, All_y_points_km, AoR_point_included, pressure_included = \
            get_points_in_aor(yaml_data, sm=sm)

        AoR.plot_aor_workflow_results(
            yaml_data, sm, All_x_points_km, All_y_points_km,
            AoR_point_included, analysis=analysis, pressure_included=pressure_included)


def set_up_aor_workflow_plots(yaml_data, well_component_type, aquifer_component_type,
                              figure_dpi=100, plot_injection_sites=False,
                              InjectionCoordx=None, InjectionCoordy=None,
                              fig_size=None, extension=None):
    """
    Sets up the plot entry dictionaries for the AoR workflow.
    """
    CriticalPressureMPa = yaml_data['Workflow']['Options'].get(
            'CriticalPressureMPa', DEFAULT_CRIT_PRESSURE_SETTING)

    brineDensityInput = yaml_data['Workflow']['Options'].get('BrineDensity', None)

    if CriticalPressureMPa != 'Calculated' and brineDensityInput is not None:
        warning_msg = ''.join([
            'The CriticalPressureMPa input provided under the Workflow: Options ',
            'section of the .yaml file was not "Calculated", but the BrineDensity ',
            'input was also provided under the same section. The BrineDensity ',
            'input is only meant to be used when CriticalPressureMPa is given as',
            ' "Calculated", so the BrineDensity input will not be used.'])

        logging.warning(warning_msg)

        # If a specific critical pressure is being used, then brineDensityInput
        # is not required. That value is only used to calculate critical pressure
        # in cases where CriticalPressureMPa == 'Calculated' and a wellbore
        # component with a brineDensity parameter is not being used. When the wellbore
        # component has a brineDensity parameter, that parameter will be used instead.
        brineDensityInput = None

    elif CriticalPressureMPa == 'Calculated' and brineDensityInput is not None:
        if well_component_type in WELL_COMPS_WITH_BRINE_DENSITY:
            warning_msg = ''.join([
                'For the AoR Workflow, the BrineDensity input (BrineDensity: {})'.format(
                    brineDensityInput),
                ' was given under the Workflow: Options section of the .yaml ',
                'file. The wellbore component used ({}) has a brine density parameter, '.format(
                    well_component_type),
                'however. The AoR Workflow will use the value of the brine ',
                'density parameter when calculating critical pressure. The ',
                'BrineDensity input provided under Workflow: Options will ',
                'be ignored. This input is meant for cases were a calculated ',
                'critical pressure is used but the wellbore component does ',
                'not have a brine density parameter.'])

            logging.warning(warning_msg)

            # Set the BrineDensity input to None. The brineDensity parameter will
            # be used instead, and removing the BrineDensity input from the .yaml
            # file set up automatically (the file that ends with '_WorkflowSetup.yaml')
            # is intended to reduce the possibility of confusion.
            brineDensityInput = None

    TimeList = None
    if 'TimeList' in yaml_data['Workflow']['Options']:
        if isinstance(yaml_data['Workflow']['Options']['TimeList'], list) or \
                yaml_data['Workflow']['Options']['TimeList'] == 'All':
            TimeList = yaml_data['Workflow']['Options']['TimeList']

        else:
            warning_msg = ''.join(['The TimeList entry was provided under ',
                                   'the Options subsection of the Workflow ',
                                   'section in the .yaml file, but the '
                                   'TimeList entry given was neither a list ',
                                   'nor the word "All". Therefore, this input ',
                                   'will be ignored.'])
            logging.warning(warning_msg)
            TimeList = None

    plot_types_to_add = AOR_PLOT_NAMES

    for plotNum, plot in enumerate(plot_types_to_add):

        for plotName in list(plot.keys()):

            # Add cases for NEW_WORKFLOWS here, if necessary
            if plotName == 'Aq_CO2_Impact_Plot':
                metric = plot[plotName].format(
                    AqCO2Metric=AOR_AQUIFER_COMPONENT_OUTPUT[
                        aquifer_component_type][0])

            elif plotName == 'Aq_Brine_Impact_Plot':
                metric = plot[plotName].format(
                    AqSaltMetric=AOR_AQUIFER_COMPONENT_OUTPUT[
                        aquifer_component_type][1])

            else:
                metric = plot[plotName]

            if extension:
                plotName = plotName + extension

            # SaveCSVFiles is required for each of the AoR plots when using
            # the AoR Workflow
            yaml_data['Plots'][plotName] = {'AoR': [metric], 'SaveCSVFiles': True,
                                            'PlotInjectionSites': plot_injection_sites,
                                            'FigureDPI': figure_dpi}

            if TimeList:
                yaml_data['Plots'][plotName]['TimeList'] = TimeList

            if metric == 'pressure':
                yaml_data['Plots'][plotName]['CriticalPressureMPa'] = CriticalPressureMPa

                if brineDensityInput:
                    yaml_data['Plots'][plotName]['BrineDensity'] = brineDensityInput

            if InjectionCoordx and InjectionCoordy:
                yaml_data['Plots'][plotName]['InjectionCoordx'] = InjectionCoordx
                yaml_data['Plots'][plotName]['InjectionCoordy'] = InjectionCoordy

            if fig_size:
                yaml_data['Plots'][plotName]['FigureSize'] = fig_size

    return yaml_data


def get_aor_aq_output_list(yaml_data, aquifer_component_type):
    """
    Returns the aquifer component outputs for the aquifer component type used
    for the AoR workflow.
    """
    output_list = AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type]

    return output_list


def set_up_aor_workflow(yaml_data):
    """
    Checks the yaml dictionary for the options required for the AoR workflow.
    """
    if 'CriticalPressureMPa' not in yaml_data['Workflow']['Options']:
        yaml_data['Workflow']['Options']['CriticalPressureMPa'] = DEFAULT_CRIT_PRESSURE_SETTING

    return yaml_data

def aor_crit_pressure_check(yaml_data, well_component_type, controls):
    """
    This function checks if the CriticalPressureMPa option was provided under
    the Options section of a .yaml file. If the entry was provided, it was
    set to a specific critical pressure value, and an OpenWellbore component is
    being used, this function makes the OpenWellbore component use that critical
    pressure.
    """
    if 'CriticalPressureMPa' in yaml_data['Workflow']['Options']:
        CriticalPressureMPa = yaml_data['Workflow']['Options']['CriticalPressureMPa']

        if CriticalPressureMPa != 'Calculated' and well_component_type == 'OpenWellbore':
            try:
                CriticalPressureMPa = float(CriticalPressureMPa)
            except:
                err_msg = ''.join([
                    'The input provided for the CriticalPressureMPa entry ({})'.format(
                        CriticalPressureMPa),
                    ', under the Workflow: Options: section of the .yaml file, ',
                    'could not be turned into a float value. Check the formatting ',
                    'used for the CriticalPressureMPa input.'])

                logging.error(err_msg)

                raise KeyError(err_msg) from None

            CriticalPressurePa = CriticalPressureMPa * 1.0e+6

            yaml_data['Workflow']['Options']['WellboreOptions']['Parameters'][
                'critPressure'] = CriticalPressurePa

            controls = {'critPressureApproach': True,
                        'enforceCritPressure': True}

    return yaml_data, controls
