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
    ControlFile_ex32a.yaml
    ControlFile_ex32b.yaml

Created: August 25th, 2022
Last Modified: February 10th, 2023

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
@contributor: Veronika Vasylkivska (Veronika.Vasylkivska@NETL.DOE.GOV)
LRST (Battelle|Leidos) supporting NETL
"""

import os
import sys
import logging
import numpy as np
import csv
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
               'pH_volume': 'pH plume volume (m$^3$)',
               'TDS_volume': 'TDS plume volume (m$^3$)',
               'Dissolved_CO2_volume': 'CO$_2$ plume volume (m$^3$)',
               'Dissolved_salt_volume': 'Salt plume volume (m$^3$)'}

TITLE_OPTIONS_LHS = {
    'pressure': ''.join(['Maximum Increase in Reservoir Pressure{}\nAcross ',
                         '{} LHS Simulations (gray: 0 MPa)']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{}\nAcross ',
                              '{} LHS Simulations (gray: 0)']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume\nAcross {} ',
                          'LHS Simulations (gray: 0 m$^3$)']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume\nAcross {} ',
                           'LHS Simulations (gray: 0 m$^3$)']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume ',
                                     'Volume\nAcross {} LHS Simulations (gray: 0 m$^3$)']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume\nAcross ',
                                      '{} LHS Simulations (gray: 0 m$^3$)'])}

TITLE_OPTIONS_FORWARD = {
    'pressure': ''.join(['Maximum Increase in Reservoir Pressure{}\nOver Time ',
                         'in the Simulation (gray: 0 MPa)']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{}\nOver ',
                              'Time in the Simulation (gray: 0)']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume\nOver Time ',
                          'in the Simulation (gray: 0 m$^3$)']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume\nOver Time ',
                           'in the Simulation (gray: 0 m$^3$)']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume Volume\nOver ',
                                     'Time in the Simulation (gray: 0 m$^3$)']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume\nOver ',
                                      'Time in the Simulation (gray: 0 m$^3$)'])}

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


def area_of_review_plot(yaml_data, model_data, output_names, sm, s,
                        output_list, locations, name='AoR_Figure1', analysis='forward',
                        savefig=None, title=None, figsize=(7, 6), figdpi=300,
                        fontname='Arial', gen_font_size=12, axis_label_font_size=14,
                        title_font_size=14, colormap=None, bold_labels=True,
                        save_results=False):
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
        Default value is (7, 6).
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
    yaml_input = get_AoR_yaml_input(yaml_data, name)

    if yaml_input['SaveCSVFiles'] is not None:
        save_results = yaml_input['SaveCSVFiles']

    if yaml_input['dpi_input'] is not None:
        figdpi = yaml_input['dpi_input']

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

        InjectionCoordx = np.array(InjectionCoordx)
        InjectionCoordx = np.unique(InjectionCoordx).tolist()
        InjectionCoordy = np.array(InjectionCoordy)
        InjectionCoordy = np.unique(InjectionCoordy).tolist()

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
                else:
                    aq_number = None
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

    if bold_labels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

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

                        if output_nm == 'pressure':
                            # Get the maximum increase in pressure in MPa
                            results[loc_ref] = max(p - values[0] for p in values) / 1e+6
                        else:
                            results[loc_ref] = max(values)

                elif analysis in ['lhs', 'parstudy']:
                    ind_list = list(range(len(time)))
                    obs_names = [full_obs_nm + '_{0}'.format(indd) for indd in ind_list]
                    obs_percentiles = percentile(s.recarray[obs_names],
                                                 [0, 25, 50, 75, 100])

                    if aq_comp_check or res_comp_check:
                        # Find the location index
                        loc_ref = int(output_component.name.split('_')[-1])

                        if output_nm == 'pressure':
                            # Get the maximum increase in pressure in MPa
                            results[loc_ref] = max(obs_percentiles[4, :]) / 1e+6
                        else:
                            results[loc_ref] = max(obs_percentiles[4, :])

    # End of loop

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
             color='k', markeredgewidth=1, markersize=6, markerfacecolor='none',
             label='Hypothetical Open Wellbore for Area of Review')
    # This is the number of columns in the legend
    ncol_number = 1

    cmap = plt.cm.get_cmap(colormap)

    if np.max(results[:, 0]) != 0:
        # I use np.max(results[:, 0]) * 1.01 here because having the maximum level
        # equal to the maximum value will cause the point with the maximum
        # value to be left out (i.e., there would be an uncolored area there).
        if np.Inf not in results[:, 0]:
            levels = np.arange(np.min(results[:, 0][results[:,0] > 0]),
                               np.max(results[:, 0]) * 1.01,
                               ((np.max(results[:, 0]) * 1.01)
                                - np.min(results[:, 0][results[:,0] > 0])) / 100)
        elif np.Inf in results[:, 0]:
            # If an OpenWellbore is placed on the injection site, you can get
            # an infinite pressure value. I include "< np.Inf" to avoid that case.
            warning_msg = ''.join([
                'The results used for the AoR plot included an infinite value. ',
                'Infinite pressures can occur if an OpenWellbore is placed ',
                'on the injection site itself. Infinite values will be ',
                'excluded from the AoR plot.'])
            logging.warning(warning_msg)
            levels = np.arange(np.min(results[:, 0][results[:,0] > 0]),
                               np.max(results[:, 0][results[:,0] < np.Inf]) * 1.01,
                               ((np.max(results[:, 0][results[:,0] < np.Inf]) * 1.01)
                                - np.min(results[:, 0][results[:,0] > 0])) / 100)

        plt.tricontourf(
            x_loc[results[:,0] < np.Inf] / 1000.0,
            y_loc[results[:,0] < np.Inf] / 1000.0,
            results[:, 0][results[:,0] < np.Inf], levels,
            locator=ticker.MaxNLocator(nbins=100, prune='lower'), cmap=colormap)

        # Plot colors for individual points so there is less ambiguity
        lgnd_check = False
        for loc_ref in range(0, len(results[:, 0])):
            if results[loc_ref, 0] != 0:
                if not lgnd_check:
                    rgba = cmap(0)
                    # This does not have any color, it's just used for the legend
                    plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000,
                              marker='o', markerfacecolor=rgba[0:3],
                              markeredgecolor='k', markeredgewidth=1.5,
                              markersize=9, linestyle='none',
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
                         markersize=9, linestyle='none')

        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))

    else:
        # If there are no results, make the colorbar go from 0 to 1
        cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                            values=np.arange(0, 1.01, 0.01),
                            format=ticker.FuncFormatter(fmt))

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
                             markeredgewidth=2, markersize=8,
                             markerfacecolor='none', zorder=1e6,
                             label='Injection Site')
                    ncol_number += 1
                else:
                    plt.plot(xcoord_val / 1000, ycoord_val / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=8,
                             markerfacecolor='none', zorder=1e6)
        else:
            plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                     marker='s', color='k', linestyle = 'none',
                     markeredgewidth=2, markersize=8,
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

    ax.legend(fancybox=False, fontsize=gen_font_size - 4, ncol=1, # ncol_number,
              edgecolor=[0, 0, 0], loc='upper center', bbox_to_anchor=(0.5, -0.25),
              framealpha=0.67)

    if output_nm in ['pressure', 'CO2saturation']:
        aq_number = ''

    # Get title corresponding to the output name
    if analysis == 'forward':
        fig_title = TITLE_OPTIONS_FORWARD.get(
            output_nm, '{} at each point (gray: zero)'.format(output_nm))

        if output_nm in TITLE_OPTIONS_FORWARD:
            fig_title = fig_title.format(aq_number)

    elif analysis in ['lhs', 'parstudy']:
        realization_number = yaml_data['ModelParams']['Analysis']['siz']
        fig_title = TITLE_OPTIONS_LHS.get(
            output_nm, 'Maximum {} at each point\nacross {} LHS simulations (gray: zero)'.format(
                output_nm, realization_number))

        if output_nm in TITLE_OPTIONS_LHS:
            fig_title = fig_title.format(aq_number, realization_number)

    # Set figure title
    plt.title(fig_title, fontsize=title_font_size,
              fontweight=selected_labelfontweight)

    plt.tight_layout()

    if title:
        plt.suptitle(title, fontweight=selected_labelfontweight,
                     fontsize=title_font_size)

    if savefig:
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

            output_dir = model_data['OutputDirectory']
            if not os.path.exists(os.path.join(output_dir, 'csv_files')):
                os.mkdir(os.path.join(output_dir, 'csv_files'))

            file_name_addition = CSV_FILE_NAME_TAGS.get(output_nm, output_nm)

            filename = os.path.join(output_dir, 'csv_files',
                                    'AoR_{}.csv'.format(file_name_addition))

            # Save the ouput for the simulation
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                for row_ref in range(0, len(x_loc) + 1):
                    writer.writerow(results_formatted[row_ref, :])

        try:
            plt.savefig(savefig, bbox_inches='tight', dpi=figdpi)
        except ValueError:
            # User has specified plot with a '.' in name but no extension.
            # Add .png as output format.
            savefig += '.png'
            plt.savefig(savefig, bbox_inches='tight', dpi=figdpi)

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
        'InjectionCoordy', 'SaveCSVFiles']

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
