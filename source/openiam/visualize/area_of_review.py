# -*- coding: utf-8 -*-
"""
Code to create a map-view figure showing an area of review (AoR) based on
metrics like pH_volume and TDS_volume. This approach is based on the work of
Bacon et al. (2020), "Probabilistic risk-based Area of Review (AoR) determination
for a deep-saline carbon storage site."

Examples illustrating applications or setup of area_of_review_plot method:
    ControlFile_ex31a.yaml
    ControlFile_ex31b.yaml
    ControlFile_ex31c.yaml
    ControlFile_ex31d.yaml
    ControlFile_ex32a.yaml
    ControlFile_ex32b.yaml
    ControlFile_ex32c.yaml

Created: August 25th, 2022
Last Modified: February 10th, 2023

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
@contributor: Veronika Vasylkivska (Veronika.Vasylkivska@NETL.DOE.GOV)
LRST (Battelle|Leidos) supporting NETL
"""
import csv
import os
import sys
import logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import cm
from matk.sampleset import percentile

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)

try:
    import openiam as iam
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


AOR_RESERVOIR_COMPONENTS = ['LookupTableReservoir', 'SimpleReservoir',
                        'AnalyticalReservoir', 'GenericReservoir']

BACKGROUND_COLOR = [0.67, 0.67, 0.67]

COLORMAP_OPTIONS = {'pressure': 'viridis',
                    'CO2saturation': 'plasma',
                    'pH_volume': 'YlGnBu',
                    'TDS_volume': 'YlOrRd',
                    'Dissolved_CO2_volume': 'YlGnBu',
                    'Dissolved_salt_volume': 'YlOrRd'}

DEFAULT_COLORMAP = 'cividis'

RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}

CBAR_LABELS = {'pressure': 'Increase in Pressure (MPa)',
               'CO2saturation': 'CO$_2$ Saturation [-]',
               'pH_volume': 'pH Plume Volume (m$^3$)',
               'TDS_volume': 'TDS Plume Volume (m$^3$)',
               'Dissolved_CO2_volume': 'CO$_2$ Plume Volume (m$^3$)',
               'Dissolved_salt_volume': 'Salt Plume Volume (m$^3$)'}

TITLE_OPTIONS_LHS = {
    'pressure': ''.join(['Maximum Increase in Reservoir Pressure{} Over Time ',
                         'Across {} LHS\nSimulations{}{}(Gray: 0 MPa)']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} Over Time ',
                              'Across {} LHS\nSimulations{}{}(Gray: 0)']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume Over Time ',
                          'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$)']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume Over Time ',
                           'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$)']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume ',
                                     'Volume Over Time Across {} LHS\n',
                                     'Simulations{}{}(Gray: 0 m$^3$)']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume ',
                                      'Over Time Across {} LHS\nSimulations',
                                      '{}{}(Gray: 0 m$^3$)'])}

TITLE_OPTIONS_FORWARD = {
    'pressure': ''.join(['Maximum Increase in Reservoir Pressure{} Over Time{}\n',
                         '{}(Gray: 0 MPa)']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} Over Time{}\n',
                              '{}(Gray: 0)']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume Over Time{}\n',
                          '{}(Gray: 0 m$^3$)']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume Over Time{}\n',
                           '{}(Gray: 0 m$^3$)']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume Volume Over ',
                                     'Time{}\n{}(Gray: 0 m$^3$)']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume Over ',
                                      'Time{}\n{}(Gray: 0 m$^3$)'])}

TITLE_OPTIONS_LHS_T_INDEX = {
    'pressure': ''.join(['Maximum Increase in Reservoir Pressure{} at t = {} years ',
                         'Across {} LHS\nSimulations{}{}(Gray: 0 MPa)']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} at t = {} years ',
                              'Across {} LHS\nSimulations{}{}(Gray: 0)']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume at t = {} years ',
                          'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$)']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume at t = {} years ',
                           'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$)']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume Volume at ',
                                     't = {} years Across {} LHS\n',
                                     'Simulations{}{}(Gray: 0 m$^3$)']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume at ',
                                      't = {} years Across {} LHS\n',
                                      'Simulations{}{}(Gray: 0 m$^3$)'])}

TITLE_OPTIONS_FORWARD_T_INDEX = {
    'pressure': ''.join(['Maximum Increase in Reservoir Pressure{} at ',
                         't = {} years{}\n{}(Gray: 0 MPa)']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} at ',
                              't = {} years{}\n{}(Gray: 0)']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume at t = {} years{}\n',
                          '{}(Gray: 0 m$^3$)']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume at t = {} years{}\n',
                           '{}(Gray: 0 m$^3$)']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume Volume ',
                                     'at t = {} years{}\n{}(Gray: 0 m$^3$)']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume ',
                                      'at t = {} years{}\n{}(Gray: 0 m$^3$)'])}

TITLE_RANGE = {
    'pressure': 'Range: {} MPa to {} MPa ',
    'CO2saturation': 'Range: {} to {} ',
    'pH_volume': 'Range: {} m$^3$ to {} m$^3$ ',
    'TDS_volume': 'Range: {} m$^3$ to {} m$^3$ ',
    'Dissolved_CO2_volume': 'Range: {} m$^3$ to {} m$^3$ ',
    'Dissolved_salt_volume': 'Range: {} m$^3$ to {} m$^3$ ',}

TITLE_SINGLE_VALUE = {
    'pressure': 'Value: {} MPa ',
    'CO2saturation': 'Value: {} ',
    'pH_volume': 'Value: {} m$^3$ ',
    'TDS_volume': 'Value: {} m$^3$ ',
    'Dissolved_CO2_volume': 'Value: {} m$^3$ ',
    'Dissolved_salt_volume': 'Value: {} m$^3$ ',}

CSV_FILE_COLUMNS = {'pressure': 'Max increase in pressure [MPa]',
                    'CO2saturation': 'Max CO2 saturation [-]',
                    'pH_volume': 'Max pH plume volume [m^3]',
                    'TDS_volume': 'Max TDS plume volume [m^3]',
                    'Dissolved_CO2_volume': 'Max CO2 plume volume [m^3]',
                    'Dissolved_salt_volume': 'Max salt plume volume [m^3]'}

CSV_FILE_NAME_TAGS = {'pressure': 'Pressure', 'CO2saturation': 'CO2_Saturation',
                      'pH_volume': 'pH_Volume', 'TDS_volume': 'TDS_Volume',
                      'Dissolved_CO2_volume': 'Dissolved_CO2_Volume',
                      'Dissolved_salt_volume': 'Dissolved_Salt_Volume'}

MAX_VAL_ADJUST = 1.1

def area_of_review_plot(yaml_data, model_data, output_names, sm, s,
                        output_list, locations, name='AoR_Figure1', analysis='forward',
                        savefig=None, title=None, figsize=(10, 8), figdpi=100,
                        fontname='Arial', gen_font_size=12, axis_label_font_size=14,
                        title_font_size=14, colormap=None, bold_labels=True,
                        save_results=False, enforce_levels=True):
    """
    Makes a map-view figure showing the maximum values of a given metric (e.g.,
    pressure, CO2saturation, pH_volume, or TDS_volume) at each location used for
    for an OpenWellbore and aquifer component (e.g., FutureGen2Aquifer). Some of
    these plots (e.g., pH_volume and TDS_volume) are meant to help define an Area
    of Review (AoR).

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param model_data: Input from the 'ModelParams' section of the .yaml file
    :type model_data: dict

    :param output_names: List of observation names to match with component models and plot.
    :type output_names: list

    :param locations: dictionary of locations assigned to each applicable component
    :type locations: dict

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param s: SampleSet with results of analysis performed on the system model.
        None for forward analysis.
    :type s: matk.SampleSet object

    :param output_list: Dictionary mapping component models to observations
    :type output_list: dict

    :param name: Figure Name to be used/created.
    :type name: str

    :param analysis: Type of OpenIAM system analysis performed ('forward',
        'lhs', or 'parstudy')
    :type analysis: str

    :param savefig: Filename to save resulting figure to. No file saved if None.
    :type savefig: str

    :param title: Optional Title for figure
    :type title: str

    :param figsize: width and height of the figure (width, height), in inches.
        Default value is (10, 8).
    :type figsize: tuple

    :param figdpi: dpi (dots-per-inch) for the figure
    :type figdpi: float or int

    :param fontname: name of the font type to be used
    :type fontname: str

    :param gen_font_size: fontsize for tick labels, etc.
    :type gen_font_size: float or int

    :param axis_label_font_size: fontsize for x and y axis labels
    :type axis_label_font_size: float or int

    :param title_font_size: fontsize for the title
    :type title_font_size: float or int

    :param colormap: string designation for a particular colormap (e.g., 'viridis')
    :type bold_labels: str

    :param bold_labels: option to use bold x and y labels and bold titles. Set to
        True for bold labels, False for normal labels.
    :type bold_labels: bool

    :param save_results: option to save the maximum values at each x and y value
        as a .csv file
    :type save_results: bool

    :return: None
    """
    # When a hypothetical wellbore is placed on the injection site itself, the
    # pressure results can include nan or Inf values. The code is set up to
    # ignore those results, but they can still cause numpy to print warning
    # statements. This statement supresses those warning statements.
    np.seterr(invalid='ignore')

    time = sm.time_array / 365.25

    yaml_input = get_AoR_yaml_input(yaml_data, name)

    if yaml_input['SaveCSVFiles'] is not None:
        save_results = yaml_input['SaveCSVFiles']

    if yaml_input['dpi_input'] is not None:
        figdpi = yaml_input['dpi_input']

    # This option specifies whether to evaluate the max. values over all times
    # (False) or evaluate the max. values for specific times (True).
    time_option = False
    time_index_list = None
    if yaml_input['TimeList'] is not None:
        time_option = True
        time_list = yaml_input['TimeList']

        if time_list == 'All':
            time_index_list = range(0, len(time))
        else:
            time_index_list = get_t_indices(time_list, time)

    components = list(sm.component_models.values())
    # If injection sites need to be plotted, get the injection sites
    if yaml_input['plot_injection_sites']:
        InjectionCoordx = []
        InjectionCoordy = []

        for comp in components:
            if comp.class_type in AOR_RESERVOIR_COMPONENTS:
                if comp.class_type != 'LookupTableReservoir':
                    InjectionCoordx.append(comp.injX)
                    InjectionCoordy.append(comp.injY)

        InjectionCoordx = np.unique(InjectionCoordx).tolist()
        InjectionCoordy = np.unique(InjectionCoordy).tolist()

    else:
        InjectionCoordx = None
        InjectionCoordy = None

    aq_number = None

    # Get the OpenWellbore data required for process_wellbore_locations()
    for comp_model in model_data['Components']:
        comp_data = yaml_data[comp_model]

        ow_cmpnt_found = False
        if 'Type' in comp_data:
            if comp_data['Type'] == 'OpenWellbore':
                comp_model_ow = comp_model
                comp_data_ow = comp_data

                if 'LeakTo' in comp_data:
                    leakTo = comp_data['LeakTo']

                    if leakTo[0:7] == 'aquifer':
                        aq_number = int(leakTo[7:None])

                ow_cmpnt_found = True
                grid_option = 'grid' in comp_data_ow['Locations']

                break

    if not ow_cmpnt_found:
        err_msg = "".join(["'Type: OpenWellbore' was not found in the components ",
                           "specified in the .yaml file. This component type is ",
                           "required for the creation of an AoR plot."])
        logging.error(err_msg)
        raise KeyError(err_msg)

    # Get locations associated with given open wellbore component
    locations_ow = locations[comp_model_ow]

    x_loc = np.array(locations_ow['coordx'])
    y_loc = np.array(locations_ow['coordy'])

    if not time_option:
        results = get_AoR_results(x_loc, output_names, sm, s, output_list,
                                  analysis=analysis, time_option=time_option)

        plot_AoR_results(aq_number, x_loc, y_loc, results, yaml_data,
                         model_data, output_names, sm, name=name, analysis=analysis,
                         savefig=savefig, title=title, figsize=figsize, figdpi=figdpi,
                         gen_font_size=gen_font_size,
                         axis_label_font_size=axis_label_font_size,
                         title_font_size=title_font_size, colormap=colormap,
                         bold_labels=bold_labels, save_results=save_results,
                         InjectionCoordx=InjectionCoordx, InjectionCoordy=InjectionCoordy,
                         grid_option=grid_option)
    else:
        min_value = None
        max_value = None

        if enforce_levels:
            min_value = 9.9e+99
            max_value = -9.9e+99

            for time_index in time_index_list:
                results = get_AoR_results(x_loc, output_names, sm, s, output_list,
                                          analysis=analysis, time_option=time_option,
                                          time_index=time_index)
                if len(results[results > 0].tolist()) > 0:
                    if np.min(results[results > 0]) < min_value:
                        min_value = np.min(results[results > 0])

                if len(results[results > 0].tolist()) > 0:
                    if np.max(results[results > 0]) > max_value:
                        results_temp = np.max(results[results > 0])
                        max_value = np.min(results_temp[results_temp < np.Inf])

            if min_value == 9.9e+99:
                min_value = 0

            if max_value == -9.9e+99:
                max_value = 1

        for time_index in time_index_list:
            results = get_AoR_results(x_loc, output_names, sm, s, output_list,
                                      analysis=analysis, time_option=time_option,
                                      time_index=time_index)

            plot_AoR_results(aq_number, x_loc, y_loc, results, yaml_data,
                             model_data, output_names, sm, name=name,
                             analysis=analysis, savefig=savefig, title=title,
                             figsize=figsize, figdpi=figdpi,
                             gen_font_size=gen_font_size,
                             axis_label_font_size=axis_label_font_size,
                             title_font_size=title_font_size, colormap=colormap,
                             bold_labels=bold_labels, save_results=save_results,
                             time_option=time_option, time_index=time_index,
                             InjectionCoordx=InjectionCoordx, InjectionCoordy=InjectionCoordy,
                             grid_option=grid_option, enforce_levels=enforce_levels,
                             min_value=min_value, max_value=max_value)


def get_AoR_results(x_loc, output_names, sm, s, output_list, analysis='forward',
                    time_option=False, time_index=None):
    """
    Evaluates and returns the maximum values of a metric for all locations.
    These maximum values are then used in a plot that is meant to inform the
    boundaries of an area of review (AoR). If time_option is set to True, the
    results will only be evaluated at the time_index provided. Otherwise, the
    results returned are the maximum values across all times.
    """
    time = sm.time_array / 365.25

    # This is used to store the maximum value of a metric at each location
    results = np.zeros((len(x_loc), 1))

    for output_nm in output_names:
        # output_list is all the components with augmented names
        for output_component in list(output_list.keys()):

            if output_nm in output_list[output_component]:
                full_obs_nm = '.'.join([output_component.name, output_nm])
                aq_comp_check = False
                res_comp_check = False

                if isinstance(output_component, (iam.FutureGen2Aquifer,
                                                 iam.FutureGen2AZMI,
                                                 iam.GenericAquifer,
                                                 iam.DeepAlluviumAquifer)):
                    aq_comp_check = True

                if isinstance(output_component, (iam.SimpleReservoir,
                                                 iam.AnalyticalReservoir,
                                                 iam.LookupTableReservoir)):
                    res_comp_check = True

                # Get maximum values at each location
                if analysis == 'forward':
                    values = sm.collect_observations_as_time_series(
                        output_component, output_nm)

                    if aq_comp_check or res_comp_check:
                        # Find the location index
                        loc_ref = int(output_component.name.split('_')[-1])

                        if not time_option:
                            if output_nm == 'pressure':
                                # Get the maximum increase in pressure in MPa.
                                results[loc_ref] = max(
                                    p - values[0] for p in values) / 1e+6
                            else:
                                results[loc_ref] = max(values)
                        else:
                            if output_nm == 'pressure':
                                # Get the maximum increase in pressure in MPa
                                results[loc_ref] = (values[time_index]
                                                    - values[0]) / 1e+6
                            else:
                                results[loc_ref] = values[time_index]

                elif analysis in ['lhs', 'parstudy']:
                    ind_list = list(range(len(time)))

                    obs_names = [full_obs_nm + '_{0}'.format(indd)
                                 for indd in ind_list]
                    obs_percentiles = percentile(s.recarray[obs_names],
                                                 [0, 25, 50, 75, 100])

                    obs_t0 = [full_obs_nm + '_0']
                    obs_percentiles_t0 = percentile(s.recarray[obs_t0],
                                                 [0, 25, 50, 75, 100])

                    if aq_comp_check or res_comp_check:
                        # Find the location index
                        loc_ref = int(output_component.name.split('_')[-1])

                        if not time_option:
                            if output_nm == 'pressure':
                                # Get the maximum increase in pressure in MPa
                                results[loc_ref] = max(
                                    obs_percentiles[4, :]
                                    - obs_percentiles_t0[4]) / 1e+6
                            else:
                                results[loc_ref] = max(obs_percentiles[4, :])
                        else:
                            if output_nm == 'pressure':
                                # Get the maximum increase in pressure in MPa
                                results[loc_ref] = (
                                    obs_percentiles[4, time_index]
                                    - obs_percentiles_t0[4]) / 1e+6
                            else:
                                results[loc_ref] = obs_percentiles[4, time_index]

    return results


def plot_AoR_results(aq_number, x_loc, y_loc, results, yaml_data, model_data,
                     output_names, sm, name='AoR_Figure1',
                     analysis='forward',savefig=None, title=None, figsize=(10, 8),
                     figdpi=100, gen_font_size=12,
                     axis_label_font_size=14, title_font_size=14, colormap=None,
                     bold_labels=True, save_results=False, time_option=False,
                     time_index=None, InjectionCoordx=None, InjectionCoordy=None,
                     grid_option=True, enforce_levels=False, min_value=None,
                     max_value=None):
    """
    Plots maximum results across all x and y values (x_loc and y_loc) for either
    all times (time_option is False) or a specific time (time option is True).
    """
    time = sm.time_array / 365.25

    if bold_labels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    yaml_input = get_AoR_yaml_input(yaml_data, name)

    output_nm = output_names[0]
    # If no colormap is specified, select a colormap. This selection
    # is done so that the default option of colormap = None produces
    # distinct colormaps for different metrics (e.g., pH_volume vs. TDS_volume).
    if colormap is None:
        colormap = COLORMAP_OPTIONS.get(output_nm, DEFAULT_COLORMAP)

    # Figures
    font = RC_FONT
    font['size'] = gen_font_size
    plt.rc('font', **font)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    # AoR Figure
    fig = plt.figure(figsize=figsize, dpi=figdpi)

    ax = plt.gca()
    ax.set_facecolor(BACKGROUND_COLOR)
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=5, markerfacecolor='none',
             label='Hypothetical Open Wellbore for Area of Review')

    # This is the number of columns in the legend
    ncol_number = 1

    cmap = plt.cm.get_cmap(colormap)

    if np.max(results[:, 0]) != 0:
        # I use np.max(results[:, 0]) * MAX_VAL_ADJUST in levels  because having
        # the maximum level equal to the maximum value can cause the area with
        # the highest values to be left out (i.e., there would be an uncolored
        # area there).
        if np.Inf not in results[:, 0]:
            min_val = np.min(results[:, 0][results[:, 0] > 0])
            a, b = '{:.2e}'.format(min_val).split('e')
            b = int(b)
            min_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            # Having > 0 here prevents the inclusion of nan values
            max_val = np.max(results[:, 0][results[:, 0] > 0])
            a, b = '{:.2e}'.format(max_val).split('e')
            b = int(b)
            max_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            if min_val == max_val:
                range_str = TITLE_SINGLE_VALUE.get(
                    output_nm, 'Value: {} ')
                range_str = range_str.format(min_val_str)

            else:
                range_str = TITLE_RANGE.get(
                    output_nm, 'Range: {} to {} ')
                range_str = range_str.format(min_val_str, max_val_str)

            if enforce_levels and time_option:
                min_val = min_value
                max_val = max_value

            interval = ((max_val * MAX_VAL_ADJUST) - min_val) / 100
            levels = np.arange(min_val, (max_val * MAX_VAL_ADJUST) + interval,
                               interval)

        elif np.Inf in results[:, 0]:
            # If an OpenWellbore is placed on the injection site, you can get
            # an infinite pressure value. I include "< np.Inf" to avoid that case.
            warning_msg = ''.join([
                'The results used for the AoR plot included an infinite value. ',
                'Infinite pressures can occur if an OpenWellbore is placed ',
                'at the injection location itself. Infinite values will be ',
                'excluded from the AoR plot.'])
            logging.warning(warning_msg)

            min_val = np.min(results[:, 0][results[:, 0] > 0])
            a, b = '{:.2e}'.format(min_val).split('e')
            b = int(b)
            min_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            max_val = np.max(results[:, 0][results[:, 0] < np.Inf])
            a, b = '{:.2e}'.format(max_val).split('e')
            b = int(b)
            max_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            if min_val == max_val:
                range_str = TITLE_SINGLE_VALUE.get(
                    output_nm, 'Value: {} ')
                range_str = range_str.format(min_val_str)

            else:
                range_str = TITLE_RANGE.get(
                    output_nm, 'Range: {} to {} ')
                range_str = range_str.format(min_val_str, max_val_str)

            if enforce_levels and time_option:
                min_val = min_value
                max_val = max_value

            interval = ((max_val * MAX_VAL_ADJUST) - min_val) / 100
            levels = np.arange(min_val, (max_val * MAX_VAL_ADJUST) + interval,
                               interval)

        try:
            # Gets rid of any nan values (from a wellbore being on the injection site)
            x_loc_temp = x_loc[results[:, 0] > 0]
            y_loc_temp = y_loc[results[:, 0] > 0]
            results_temp = results[results[:, 0] > 0]

            # Gets rid of any Inf values (from a wellbore being on the injection site)
            x_loc_temp = x_loc[results[:, 0] < np.Inf]
            y_loc_temp = y_loc[results[:, 0] < np.Inf]
            results_temp = results[results[:, 0] < np.Inf]

            plt.tricontourf(
                x_loc_temp / 1000.0, y_loc_temp / 1000.0,
                results_temp[:, 0], levels, cmap=colormap,
                locator=ticker.MaxNLocator(nbins=100, prune='lower'))

            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

        except:
            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

        # Plot colors for individual points so there is less ambiguity
        lgnd_check = False
        for loc_ref in range(0, len(results[:, 0])):
            if results[loc_ref, 0] > 0 and results[loc_ref, 0] < np.Inf:
                if not lgnd_check:
                    rgba = cmap(0)
                    # This does not have any color, it's just used for the legend
                    plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000,
                              marker='o', markerfacecolor=rgba[0:3],
                              markeredgecolor='k', markeredgewidth=1.5,
                              markersize=12, linestyle='none',
                              label='Wellbore with Nonzero Result')
                    lgnd_check = True
                    ncol_number += 1

                # Get the color by the value's proximity to the upper and lower
                # limits used for levels (which set the colorbar limits above)
                rgba = cmap((results[loc_ref, 0] - np.min(levels))
                             / (np.max(levels) - np.min(levels)))
                plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000,
                         marker='o', markerfacecolor=rgba[0:3],
                         markeredgecolor='k', markeredgewidth=1.5,
                         markersize=12, linestyle='none')

    else:
        if enforce_levels and time_option:
            min_val = min_value
            max_val = max_value

            interval = ((max_val * MAX_VAL_ADJUST) - min_val) / 100
            levels = np.arange(min_val, (max_val * MAX_VAL_ADJUST) + interval,
                               interval)

            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

            range_str = ''

        else:
            # If there are no results and no enforced colorbar, make the
            # colorbar go from 0 to 1.
            levels = np.arange(0, 1.01, 0.01)
            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

            range_str = ''

    # Default cbar label is output name
    cbar_label = CBAR_LABELS.get(output_nm, output_nm)
    cbar.set_label(
        cbar_label, rotation=90, fontsize=axis_label_font_size,
        fontweight=selected_labelfontweight)

    if yaml_input['plot_injection_sites']:
        if len(InjectionCoordx) == 0 and yaml_input['InjectionCoordx'] is not None:
            # The lists will be empty for LookupTableReservoirs, use the .yaml input
            InjectionCoordx = yaml_input['InjectionCoordx']
            InjectionCoordy = yaml_input['InjectionCoordy']

        if isinstance(InjectionCoordx, list) and len(InjectionCoordx) > 0:
            for injRef, (xcoord_val, ycoord_val) in enumerate(
                    zip(InjectionCoordx, InjectionCoordy)):
                if injRef == 0:
                    plt.plot(xcoord_val / 1000, ycoord_val / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=6,
                             markerfacecolor='none', zorder=1e6,
                             label='Injection Site')
                    ncol_number += 1
                else:
                    plt.plot(xcoord_val / 1000, ycoord_val / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=6,
                             markerfacecolor='none', zorder=1e6)
        else:
            plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                     marker='s', color='k', linestyle = 'none',
                     markeredgewidth=2, markersize=6,
                     markerfacecolor='none', zorder=1e6,
                     label='Injection Site')
            ncol_number += 1

    x_loc_copy = x_loc.tolist()
    y_loc_copy = y_loc.tolist()
    if yaml_input['plot_injection_sites']:
        x_loc_copy = x_loc.tolist()
        y_loc_copy = y_loc.tolist()
        if isinstance(InjectionCoordx, list) and len(InjectionCoordx) > 0:
            for (xcoord_val, ycoord_val) in zip(InjectionCoordx, InjectionCoordy):
                x_loc_copy.append(xcoord_val)
                y_loc_copy.append(ycoord_val)
        else:
            x_loc_copy.append(InjectionCoordx)
            y_loc_copy.append(InjectionCoordy)

    # These contain the OpenWellbore locations and the injection locations,
    # if the injection locations are being plotted.
    x_loc_copy = np.array(x_loc_copy)
    y_loc_copy = np.array(y_loc_copy)

    if grid_option:
        # These are used for the buffer room b/c they will not include the
        # injection location(s), only the grid-based OpenWellbore locations.
        x_vals_temp = np.unique(x_loc)
        y_vals_temp = np.unique(y_loc)
        plt.xlim((np.min(x_loc_copy) - ((x_vals_temp[1]
                                         - x_vals_temp[0]))) / 1000.0,
                  (np.max(x_loc_copy) + ((x_vals_temp[1]
                                          - x_vals_temp[0]))) / 1000.0)
        plt.ylim((np.min(y_loc_copy) - ((y_vals_temp[1]
                                         - y_vals_temp[0]))) / 1000.0,
                  (np.max(y_loc_copy) + ((y_vals_temp[1]
                                          - y_vals_temp[0]))) / 1000.0)
    else:
        xlim_adjust_val = (np.max(x_loc_copy) - np.min(x_loc_copy)) / 20
        ylim_adjust_val = (np.max(y_loc_copy) - np.min(y_loc_copy)) / 20
        plt.xlim((np.min(x_loc_copy) - xlim_adjust_val) / 1000.0,
                 (np.max(x_loc_copy) + xlim_adjust_val) / 1000.0)
        plt.ylim((np.min(y_loc_copy) - ylim_adjust_val) / 1000.0,
                 (np.max(y_loc_copy) + ylim_adjust_val) / 1000.0)

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    ax.set_aspect('equal')
    ax.set_xlabel('Easting (km)', fontsize=axis_label_font_size,
                  fontweight=selected_labelfontweight)
    ax.set_ylabel('Northing (km)', fontsize=axis_label_font_size,
                  fontweight=selected_labelfontweight)

    ax.legend(fancybox=False, fontsize=gen_font_size - 2, ncol=ncol_number,
              edgecolor=[0, 0, 0], loc='upper center', bbox_to_anchor=(0.5, -0.1),
              framealpha=0.67)

    if output_nm in ['pressure', 'CO2saturation']:
        aq_number = ''

    # Get title corresponding to the output name
    if analysis == 'forward':
        if not time_option:
            # Need a space afterwards for the case where fig_title is set to
            # '{} at Each Point{}{}(Gray: Zero)'
            if range_str != '':
                comma_str = ', '
            else:
                comma_str = ' '

            fig_title = TITLE_OPTIONS_FORWARD.get(
                output_nm, '{} at Each Point{}{}(Gray: Zero)'.format(
                    output_nm, comma_str, range_str))

            # No space afterwards b/c the comma_str goes right before a '\n'
            if range_str != '':
                comma_str = ','
            else:
                comma_str = ''

            if output_nm in TITLE_OPTIONS_FORWARD:
                fig_title = fig_title.format(aq_number, comma_str, range_str)
        else:
            # No space afterwards b/c the comma_str goes right before a '\n'
            if range_str != '':
                comma_str = ','
            else:
                comma_str = ''

            fig_title = TITLE_OPTIONS_FORWARD_T_INDEX.get(
                output_nm, '{} at Each Point at Time\nt = {} years{}{}(Gray: Zero)'.format(
                    output_nm, time[time_index], comma_str, range_str))

            if output_nm in TITLE_OPTIONS_FORWARD_T_INDEX:
                fig_title = fig_title.format(aq_number, time[time_index],
                                             comma_str, range_str)

    elif analysis in ['lhs', 'parstudy']:
        realization_number = yaml_data['ModelParams']['Analysis']['siz']

        if range_str != '':
            comma_str = ', '
        else:
            comma_str = ' '

        if not time_option:
            fig_title = TITLE_OPTIONS_LHS.get(
                output_nm,
                ''.join(['Maximum {} at Each Point\nAcross {} ',
                         'LHS Simulations{}{}(Gray: Zero)']).format(
                    output_nm, realization_number, comma_str, range_str))

            if output_nm in TITLE_OPTIONS_LHS:
                fig_title = fig_title.format(aq_number, realization_number,
                                             comma_str, range_str)
        else:
            fig_title = TITLE_OPTIONS_LHS_T_INDEX.get(
                output_nm,
                ''.join(['Maximum {} at Each Point at Time t = {} years\n',
                         'Across {} LHS Simulations{}{}(Gray: Zero)']).format(
                    output_nm, time[time_index], realization_number,
                    comma_str, range_str))

            if output_nm in TITLE_OPTIONS_LHS_T_INDEX:
                fig_title = fig_title.format(aq_number, time[time_index],
                                             realization_number,
                                             comma_str, range_str)

    # Set figure title
    plt.title(fig_title, fontsize=title_font_size,
              fontweight=selected_labelfontweight)

    plt.subplots_adjust(left=0.001, bottom=0.15, right=0.901,
                        top=0.875, wspace=0.1, hspace=0.1)
    # plt.tight_layout()

    if title:
        plt.suptitle(title, fontweight=selected_labelfontweight,
                     fontsize=title_font_size)

    if savefig:
        output_dir = model_data['OutputDirectory']

        if save_results:
            results_formatted = np.empty(((len(x_loc) + 1), 3), dtype=list)

            results_formatted[0, 0] = 'x (km)'
            results_formatted[0, 1] = 'y (km)'

            for row_ref, (x_val, y_val) in enumerate(zip(x_loc, y_loc)):
                results_formatted[row_ref + 1, 0] = str(x_val / 1000)
                results_formatted[row_ref + 1, 1] = str(y_val / 1000)

            # Get the column name depending on the output
            # If name does not match known observations use it as a column label
            results_formatted[0, 2] = CSV_FILE_COLUMNS.get(output_nm, output_nm)

            results_formatted[1:None, 2:None] = results

            if not os.path.exists(os.path.join(output_dir, 'csv_files')):
                os.mkdir(os.path.join(output_dir, 'csv_files'))

            file_name_addition = CSV_FILE_NAME_TAGS.get(output_nm, output_nm)

            if not time_option:
                filename = os.path.join(output_dir, 'csv_files',
                                        'AoR_{}.csv'.format(file_name_addition))
            else:
                filename = os.path.join(output_dir, 'csv_files',
                                        'AoR_{}_tIndex_{:.0f}.csv'.format(
                                            file_name_addition, time_index))

            # Save the ouput for the simulation
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                for row_ref in range(0, len(x_loc) + 1):
                    writer.writerow(results_formatted[row_ref, :])

            f.close()

        if '.' in name:
            name_main = name[0:name.index('.')]
            name_main += '_tIndex_{:.0f}'.format(time_index)
            name_extension = name[name.index('.'):]
        else:
            name_main = name
            name_extension = '.png'

        if time_option:
            name_main += '_tIndex_{:.0f}'.format(time_index)

        plt.savefig(os.sep.join([output_dir,
                                 name_main + name_extension]), dpi=figdpi)
        plt.close()
    else:
        plt.show()


def not_boolean_debug_message_AoR(input_name, name, default_value):
    """
    Returns string delivering debug message regarding a variable not being
    of boolean type and setting it to the default value (True or False).

    input_name: string
    default_value: True or False
    """
    msg = ''.join(['Input provided for {} within the AoR plot {} was not of ',
                   'boolean type. Using the default value of {}.']).format(
                       input_name, name, default_value)
    return msg


def get_AoR_yaml_input(yaml_data, name):
    """
    Function that reads the AoR plot input provided in the .yaml file
    and returns a dictionary containing the input. Note that InjectionCoordx
    and InjectionCoordy are only required if a LookupTableReservoir is being
    used. Other reservoir components will have locX and locY values that can
    be used to plot injection locations, but LookupTableReservoir components do
    not have locX and locY values.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param name: Name of the AoR plot provided in the Plots section of the
        .yaml file
    :type name: str

    :returns: yaml_input
    """

    yaml_input_keys = [
        'dpi_input', 'plot_injection_sites', 'InjectionCoordx',
        'InjectionCoordy', 'SaveCSVFiles', 'TimeList']

    InjectionCoord_debug_msg = ''.join([
        'InjectionCoord{} was provided for the AoR plot ', name,
        ', but not InjectionCoord{}. Check your input. Injection ',
        'sites will not be displayed.'])

    # Initialize output
    yaml_input = {key: None for key in yaml_input_keys}

    # Get shortcut to data to be analyzed
    AoR_plot_data = yaml_data['Plots'][name]

    if AoR_plot_data is not None:
        if 'FigureDPI' in AoR_plot_data:
            yaml_input['dpi_input'] = AoR_plot_data['FigureDPI']

        if 'TimeList' in AoR_plot_data:
            yaml_input['TimeList'] = AoR_plot_data['TimeList']

        if 'PlotInjectionSites' in AoR_plot_data:
            yaml_input['plot_injection_sites'] = AoR_plot_data[
                'PlotInjectionSites']
            if not isinstance(yaml_input['plot_injection_sites'], bool):
                debug_msg = not_boolean_debug_message_AoR(
                    'PlotInjectionSites', name, False)
                logging.debug(debug_msg)
                yaml_input['plot_injection_sites'] = False

        if 'InjectionCoordx' in AoR_plot_data:
            try:
                yaml_input['InjectionCoordx'] = float(AoR_plot_data['InjectionCoordx'])
            except TypeError:
                yaml_input['InjectionCoordx'] = list(AoR_plot_data['InjectionCoordx'])

        if 'InjectionCoordy' in AoR_plot_data:
            try:
                yaml_input['InjectionCoordy'] = float(AoR_plot_data['InjectionCoordy'])
            except TypeError:
                yaml_input['InjectionCoordy'] = list(AoR_plot_data['InjectionCoordy'])

            if yaml_input['InjectionCoordx'] is None:
                debug_msg = InjectionCoord_debug_msg.format('y', 'x')
                logging.debug(debug_msg)
                yaml_input['plot_injection_sites'] = False

            elif isinstance(yaml_input['InjectionCoordx'], list) and \
                    isinstance(yaml_input['InjectionCoordy'], list):
                if len(yaml_input['InjectionCoordx']) != len(
                        yaml_input['InjectionCoordy']):
                    debug_msg = ''.join([
                        'The InjectionCoordy provided for the AoR plot ', name,
                        ' was not of the same length as the InjectionCoordx ',
                        'provided. Check your input. Injection sites will not ',
                        'be displayed.'])
                    logging.debug(debug_msg)
                    yaml_input['plot_injection_sites'] = False

        if yaml_input['InjectionCoordx'] is not None and yaml_input[
                'InjectionCoordy'] is None:
            debug_msg = InjectionCoord_debug_msg.format('x', 'y')
            logging.debug(debug_msg)
            yaml_input['plot_injection_sites'] = False

        if 'SaveCSVFiles' in AoR_plot_data:
            yaml_input['SaveCSVFiles'] = AoR_plot_data['SaveCSVFiles']

            if not isinstance(yaml_input['SaveCSVFiles'], bool):
                debug_msg = not_boolean_debug_message_AoR(
                    'SaveCSVFiles', name, True)
                logging.debug(debug_msg)
                yaml_input['SaveCSVFiles'] = None

    return yaml_input


def get_t_indices(time_list, time_array):
    """
    Returns the time index corresponding with the time_array value closest to
    the times in time_list. Both time_list and time_array have units of years.
    """
    time_index_list = []

    corr_t_index = np.arange(0, len(time_array))

    if not isinstance(time_list, list):
        time_list = [time_list]

    for time in time_list:
        abs_diff = np.zeros(len(time_array))

        for t_ref, t_val in enumerate(time_array):
            abs_diff[t_ref] = np.abs(time - (t_val))

        closest_t_index = corr_t_index[abs_diff == np.min(abs_diff)]

        # If more than one time_array value had the same distance to the time
        # in question (e.g., one value is 0.5 yrs before and another is 0.5 yrs
        # after), then just pick one.
        if len(closest_t_index) > 1:
            closest_t_index = closest_t_index[-1]

        if isinstance(closest_t_index, list):
            closest_t_index = closest_t_index[0]

        time_index_list.append(int(closest_t_index))

    return time_index_list
