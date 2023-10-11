"""
Code to create figures showing the stratigraphy created with
strata.py in the folder cfi.

Examples illustrating applications or setup of stratigraphy_plot method:
    ControlFile_ex32b.yaml
    ControlFile_ex33a.yaml
    ControlFile_ex33b.yaml
    ControlFile_ex34.yaml
    ControlFile_ex35.yaml
    ControlFile_ex36.yaml
    ControlFile_ex38c.yaml

Created: September 15th, 2022
Last Modified: February 10th, 2023

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
@contributor: Veronika Vasylkivska (Veronika.Vasylkivska@NETL.DOE.GOV)
LRST (Battelle|Leidos) supporting NETL
"""

import os
import sys
import logging

import numpy as np
import matplotlib.pyplot as plt

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)

import openiam.cfi.strata as iam_strata
from openiam.visualize.get_strat_obs import get_strat_comp_obs

DEPTH_LABEL = '{:.2f} m'


def stratigraphic_column(yaml_data, sm, name='strat_column_Figure1',
                         savefig=None, title=None, figsize=(6, 10),
                         figdpi=100, genfontsize=12, axislabelfontsize=14,
                         titlefontsize=14, boldlabels=True):
    """
    Makes a figure showing the domain's stratigraphy at a particular location.
    While the figures made by stratigraphy_plot.py use a 3-dimensional perspective
    and can show a large area, the figure created by this function is a
    stratigraphic column representing one location.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param name: Figure Name to be used/created.
    :type name: str

    :param savefig: Filename to save resulting figure to. No file saved if None.
    :type savefig: str

    :param title: Optional Title for figure
    :type title: str

    :param figsize: width and height of the figure (width, height) in inches.
        Default value is (10, 6).
    :type figsize: tuple

    :param figdpi: dpi (dots-per-inch) for the figure
    :type figdpi: float or int

    :param genfontsize: fontsize for tick labels, etc.
    :type genfontsize: float or int

    :param axislabelfontsize: fontsize for x and y axis labels
    :type axislabelfontsize: float or int

    :param titlefontsize: fontsize for the title
    :type titlefontsize: float or int

    :return: None
    """
    # These lists indicate the stratigraphy component types that offer thicknesses
    # and depths as parameters or as observations.
    types_strata_pars = iam_strata.get_comp_types_strata_pars()
    types_strata_obs = iam_strata.get_comp_types_strata_obs()

    if boldlabels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    # Get the stratigrapy information from the .yaml file
    comp_type = iam_strata.get_strata_type_from_yaml(yaml_data)

    # Get the stratigraphy plot input from the .yaml file
    strat_plot_yaml_input = read_strat_col_plot_yaml_input(yaml_data, name)

    if strat_plot_yaml_input['XValue'] is not None:
        x_value = float(strat_plot_yaml_input['XValue'])
    else:
        x_value = 0

    if strat_plot_yaml_input['YValue'] is not None:
        y_value = float(strat_plot_yaml_input['YValue'])
    else:
        y_value = 0

    if strat_plot_yaml_input['dpi_input'] is not None:
        figdpi = float(strat_plot_yaml_input['dpi_input'])

    if strat_plot_yaml_input['plot_depth_text'] is not None:
        plot_depth_text = strat_plot_yaml_input['plot_depth_text']
    else:
        plot_depth_text = True

    # Get the data
    if comp_type in types_strata_pars:

        components = list(sm.component_models.values())

        for comp in components:
            if comp.class_type == comp_type:
                strata_comp = comp
                break

        numShaleLayers = strata_comp.get_num_shale_layers()

        strata_dict = iam_strata.get_strata_info_from_component(strata_comp)

        shaleThicknessList = strata_dict['shaleThicknesses']
        aquiferThicknessList = strata_dict['aquiferThicknesses']
        reservoirThickness = strata_dict['reservoirThickness']

        # These need to have negative values
        shaleTopDepthList = (-1 * np.array(strata_dict['shaleTopDepths'])).tolist()
        aquiferTopDepthList = (-1 * np.array(strata_dict['aquiferTopDepths'])).tolist()
        reservoirDepth = (-1 * np.array(strata_dict['reservoirDepth'])).tolist()
        reservoirTopDepth = (-1 * np.array(strata_dict['reservoirTopDepth'])).tolist()

    elif comp_type in types_strata_obs:
        # Get the number of shale layers
        components = list(sm.component_models.values())
        for comp in components:
            if comp.class_type in types_strata_obs:
                numShaleLayers = comp.get_num_shale_layers()

                break

        stratigraphy_by_loc_temp = get_strat_comp_obs(
            x_value, y_value, numShaleLayers, yaml_data)

        shaleThicknessList = []
        aquiferThicknessList = []
        shaleTopDepthList = []
        aquiferTopDepthList = []
        for shaleRef in range(numShaleLayers):
            shaleThicknessList.append(np.max(stratigraphy_by_loc_temp[
                'shale{}Thickness'.format(shaleRef + 1)]))

            # Must be a negative value
            shaleTopDepthList.append(-np.max(stratigraphy_by_loc_temp[
                'shale{}TopDepth'.format(shaleRef + 1)]))

            if (shaleRef + 1) < numShaleLayers:
                aquiferThicknessList.append(np.max(stratigraphy_by_loc_temp[
                    'aquifer{}Thickness'.format(shaleRef + 1)]))

                # Must be a negative value
                aquiferTopDepthList.append(-np.max(stratigraphy_by_loc_temp[
                    'aquifer{}TopDepth'.format(shaleRef + 1)]))

        reservoirThickness = np.max(stratigraphy_by_loc_temp['reservoirThickness'])

        # Must be negative values
        reservoirTopDepth = -np.max(stratigraphy_by_loc_temp['reservoirTopDepth'])
        reservoirDepth = -np.max(stratigraphy_by_loc_temp['reservoirDepth'])

    reservoirColor, reservoirAlpha, reservoirAlphaFill, reservoirLabel, \
        shaleColor, shaleAlpha, shaleAlphaFill, shaleLabel, \
            aquiferColor, aquiferAlpha, aquiferAlphaFill, aquiferLabel, \
                _, _, _, _ = \
                    iam_strata.get_plotting_setup_for_units(
                        strat_plot_yaml_input, numShaleLayers,
                        reservoirThickness=reservoirThickness,
                        shaleThicknessList=shaleThicknessList,
                        aquiferThicknessList=aquiferThicknessList,
                        include_thickness=True)


    # Make the figure
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': genfontsize}
    plt.rc('font', **font)

    # Set up of the x axis
    x_value_strat_col = 1

    # Keeping these separate in case they need to be changed in the future
    x_range_shales = [x_value_strat_col - 2, x_value_strat_col + 1.5]
    x_range_aquifers = [x_value_strat_col - 2, x_value_strat_col + 2]
    x_range_reservoir = [x_value_strat_col - 2, x_value_strat_col + 2]

    # x text locations for unit labels
    x_text_shales = x_range_shales[0] + (
        0.05 * (x_value_strat_col - x_range_shales[0]))
    x_text_aquifers = x_range_aquifers[0] + (
        0.05 * (x_value_strat_col - x_range_aquifers[0]))
    x_text_reservoir = x_range_reservoir[0] + (
        0.05 * (x_value_strat_col - x_range_reservoir[0]))

    # x text locations for unit depths
    x_text_depth_shales = x_range_shales[1] + (
        0.05 * (x_value_strat_col - x_range_shales[0]))
    x_text_depth_aquifers = x_range_aquifers[1] + (
        0.05 * (x_value_strat_col - x_range_aquifers[0]))
    x_text_depth_reservoir = x_range_reservoir[1] + (
        0.05 * (x_value_strat_col - x_range_reservoir[0]))

    if plot_depth_text:
        x_lims_strat_col = [x_value_strat_col - 3, x_value_strat_col + 3.5]
    else:
        x_lims_strat_col = [x_value_strat_col - 3, x_value_strat_col + 3]

    # Scale the buffer by the thinnest unit so that it is uniform
    y_buffer_text = np.min(aquiferThicknessList + shaleThicknessList) * 0.075

    fig = plt.figure(figsize=figsize, dpi=figdpi)
    ax = fig.add_subplot()

    ax.plot(x_range_reservoir, [reservoirDepth, reservoirDepth],
            color=reservoirColor, alpha=reservoirAlpha, zorder=10)

    ax.plot(x_range_reservoir, [reservoirTopDepth, reservoirTopDepth],
            color=reservoirColor, alpha=reservoirAlpha, zorder=10)

    ax.fill_between(x_range_reservoir, [reservoirTopDepth, reservoirTopDepth],
                    [reservoirDepth, reservoirDepth],
                    color=reservoirColor, alpha=reservoirAlphaFill)

    # Text for unit label
    y_text = reservoirDepth + y_buffer_text

    ax.text(x_text_reservoir, y_text, reservoirLabel, color=reservoirColor,
            alpha=reservoirAlpha, fontsize=genfontsize - 2,
            fontweight=selected_labelfontweight, zorder=10)

    # Text for unit depth
    if plot_depth_text:
        depth_text = DEPTH_LABEL.format(reservoirDepth)

        ax.text(x_text_depth_reservoir, y_text, depth_text,
                color=reservoirColor, alpha=reservoirAlpha, fontsize=genfontsize - 2,
                fontweight=selected_labelfontweight, zorder=10)

    for shaleRef, shaleTDValue in enumerate(shaleTopDepthList):
        ax.plot(x_range_shales, [reservoirTopDepth, reservoirTopDepth],
                color=shaleColor[shaleRef], alpha=shaleAlpha[shaleRef], zorder=11)

        ax.plot(x_range_shales, [shaleTDValue, shaleTDValue],
                color=shaleColor[shaleRef], alpha=shaleAlpha[shaleRef], zorder=11)

        if shaleRef == 0:
            ax.fill_between(
                x_range_shales, [shaleTDValue, shaleTDValue],
                [reservoirTopDepth, reservoirTopDepth],
                color=shaleColor[shaleRef], alpha=shaleAlphaFill[shaleRef])

            y_text = reservoirTopDepth + y_buffer_text
            depth_text = DEPTH_LABEL.format(reservoirTopDepth)
        else:
            ax.fill_between(
                x_range_shales, [shaleTDValue, shaleTDValue],
                [aquiferTopDepthList[shaleRef - 1],
                 aquiferTopDepthList[shaleRef - 1]],
                color=shaleColor[shaleRef], alpha=shaleAlphaFill[shaleRef])

            y_text = aquiferTopDepthList[shaleRef - 1] + y_buffer_text
            depth_text = DEPTH_LABEL.format(aquiferTopDepthList[shaleRef - 1])

        # Text for unit label
        ax.text(x_text_shales, y_text, shaleLabel[shaleRef],
                color=shaleColor[shaleRef],
                alpha=shaleAlpha[shaleRef], fontsize=genfontsize - 2,
                fontweight=selected_labelfontweight, zorder=11)

        # Text for unit depth
        if plot_depth_text:
            ax.text(x_text_depth_shales, y_text, depth_text,
                    color=shaleColor[shaleRef],
                    alpha=shaleAlpha[shaleRef], fontsize=genfontsize - 2,
                    fontweight=selected_labelfontweight, zorder=10)

        if shaleRef != (len(shaleTopDepthList) - 1):
            # Plot aquifers
            ax.plot(x_range_aquifers, [shaleTDValue, shaleTDValue],
                    color=aquiferColor[shaleRef],
                    alpha=aquiferAlpha[shaleRef], zorder=10)

            ax.plot(x_range_aquifers,
                    [aquiferTopDepthList[shaleRef], aquiferTopDepthList[shaleRef]],
                    color=aquiferColor[shaleRef], alpha=aquiferAlpha[shaleRef],
                    zorder=10)

            ax.fill_between(
                x_range_aquifers,
                [aquiferTopDepthList[shaleRef], aquiferTopDepthList[shaleRef]],
                [shaleTDValue, shaleTDValue],
                color=aquiferColor[shaleRef], alpha=aquiferAlphaFill[shaleRef])

            # Text for unit label
            y_text = shaleTopDepthList[shaleRef] + y_buffer_text

            ax.text(x_text_aquifers, y_text, aquiferLabel[shaleRef],
                    color=aquiferColor[shaleRef], alpha=aquiferAlpha[shaleRef],
                    fontsize=genfontsize - 2, fontweight=selected_labelfontweight,
                    zorder=10)

            # Text for unit depth
            if plot_depth_text:
                depth_text = DEPTH_LABEL.format(shaleTDValue)

                ax.text(x_text_depth_aquifers, y_text, depth_text,
                        color=aquiferColor[shaleRef], alpha=aquiferAlpha[shaleRef],
                        fontsize=genfontsize - 2, fontweight=selected_labelfontweight,
                        zorder=10)

    if plot_depth_text:
        ax.text(x_text_depth_shales, 0, '0 m', color=shaleColor[shaleRef],
                alpha=shaleAlpha[shaleRef], fontsize=genfontsize - 2,
                fontweight=selected_labelfontweight, zorder=10)

    ax.set_xlim(x_lims_strat_col)
    ax.set_xticks([])

    if comp_type in types_strata_pars:
        x_label = 'Stratigraphic Column'
    else:
        x_label = 'Stratigraphic Column at\nx = {} km, y = {} km'.format(
            x_value / 1000, y_value / 1000)

    ax.set_xlabel(x_label, fontsize=axislabelfontsize, fontweight=selected_labelfontweight)

    ax.set_ylabel('Depth (m)', fontsize=axislabelfontsize, fontweight='bold')

    ax.set_title('Stratigraphy for the Study Area',
                 fontsize=titlefontsize, fontweight=selected_labelfontweight)

    if title:
        fig.suptitle(title, fontweight=selected_labelfontweight,
                     fontsize=titlefontsize)

    fig.tight_layout()

    if savefig:
        try:
            fig.savefig(savefig, bbox_inches='tight', dpi=figdpi)
        except ValueError:
            # User has specified plot with a '.' in name but no extension.
            # Add .png as output format.
            savefig += '.png'
            fig.savefig(savefig, bbox_inches='tight', dpi=figdpi)
    else:
        plt.show()


def read_strat_col_plot_yaml_input(yaml_data, name):
    """
    Function that reads the Stratigraphy plot input provided in the .yaml file
    and returns a dictionary containing the input.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param name: Name of the Stratigraphy plot provided in the Plots section of
        the .yaml file
    :type name: str

    :returns: yaml_input
    """
    # Default values
    yaml_input_keys = [
        'dpi_input', 'XValue', 'YValue', 'plot_depth_text']

    # Initialize output
    yaml_input = {key: None for key in yaml_input_keys}

    # Get shortcut to data to be analyzed
    strata_plot_data = yaml_data['Plots'][name]['StratigraphicColumn']

    if strata_plot_data is not None:
        # Checks for and adds unit colors, unit alpha values, and unit labels
        # provided in the plot entry. Excludes any invalid input provided.
        yaml_input = iam_strata.check_color_alpha_label_yaml_input(
            yaml_input, strata_plot_data, name)

        if 'FigureDPI' in strata_plot_data:
            yaml_input['dpi_input'] = strata_plot_data['FigureDPI']

        if 'XValue' in strata_plot_data:
            yaml_input['XValue'] = strata_plot_data['XValue']

        if 'YValue' in strata_plot_data:
            yaml_input['YValue'] = strata_plot_data['YValue']

        if 'DepthText' in strata_plot_data:
            if isinstance(strata_plot_data['DepthText'], bool):
                yaml_input['plot_depth_text'] = strata_plot_data['DepthText']
            else:
                debug_msg = ''.join([
                    'The input provided for DepthText in the StratigraphicColumn ',
                    'plot ', name, ' was not of boolean type. The default value '
                    'of True will be used.'])
                logging.debug(debug_msg)

    return yaml_input
