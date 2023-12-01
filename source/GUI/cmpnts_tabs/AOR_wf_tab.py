"""
Module contains several methods needed for creating tab (page) in GUI
for AOR Workflow component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import BooleanVar

from dictionarydata import componentVars, componentChoices
from dictionarydata import workflowVars, workflowChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)

from cmpnts_tabs.locations import (add_inj_well_frame_widgets,
                                   add_obs_locs_frame_widgets,
                                   read_obs_locations_data)
from cmpnts_tabs.commons import commons_read_tab_vars
from cmpnts_tabs.parameter_entry import ParameterEntry

LABEL_WIDTH = 17

AOR_PARAMETERS = ['CriticalPressureMPa',
                  'BrineDensity',
                  'FigureDPI']

AOR_PARAMETERS_SETUP = {
    'CriticalPressureMPa': ["CriticalPressure [MPa]:",
                            'CriticalPressureMPa'],
    'BrineDensity': ["Brine density [kg/m{}]:".format(u'\u00B3'),
                     'Brine density'],
    'FigureDPI': ["Figure Resolution [DPI]:", "FigureDPI"]
    }

# Set AOR Workflow parameters names and value, min, max, second value, mean, std, bounds
AOR_PARAMETER_VALUES = {
    'CriticalPressureMPa': [5.0, 1.0, 20.0, 10.0, 5.0, 0.01, 0.1, 90.0],
    'BrineDensity': [1000, 900, 1500, 1200, 1200, 150, 900, 1500],
    'FigureDPI': [100, 100, 400, 100, 200, 50, 100, 400],
    }

AOR_OBSERVATIONS = []

AOR_OBSERVATIONS_SETUP = {}

brineDensityWellbores = ['Cemented Wellbore', 'Cemented Wellbore WR']


def disable_injection_location_widgets(controller, frame):
    """
    Disable injection location widgets
    on the AOR Workflow Tab based on
    whether to plot injection location.
    Only applicable in the case of the
    LookupTableReservoir.
    """

    if controller.component_list['reservoir'][1].get() == 'Lookup Table Reservoir':
        frame_state = {1: 'normal', 0: 'disabled'}
        check_button_state = frame.checkbox_variable.get()

        # Widgets for input of injection well x and y
        for field in frame.coords_fields:
            field.configure(state=frame_state[check_button_state])
            field.configure(state=frame_state[check_button_state])


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'AoR', parameter_names=AOR_PARAMETERS,
        observation_names=None)

    if 'LookupTableReservoir1' in componentVars.keys():
        cmpnt_data = read_obs_locations_data(cmpnt_data, cmpnt_nm)

    return cmpnt_data


def setup_parameter_frame(controller, frame, par_name, par_bounds,
                          par_label_text, tool_tip_text,
                          par_name_label_width, distr_menu_width,
                          distr_arg_label_width, text_field_width,
                          distr_options, par_args_variables, cmpnt_nm, tool_tip):
    """ Create, configure and place widgets for a given parameter frame."""
    # Update attributes of frame
    frame.labelText = par_label_text
    frame.component = cmpnt_nm
    frame.distType = par_args_variables['distribution']
    frame.paramVars = par_args_variables
    frame.toolTipText = tool_tip_text
    frame.text = par_name
    frame.par_bounds = {'lower_bound': None,
                        'upper_bound': None,
                        'discrete_bounds': None}
    frame.distr_options = distr_options

    # Copy provided values
    for key in par_bounds:
        frame.par_bounds[key] = par_bounds[key]

    # Create common widgets
    # Label with parameter name and units
    par_name_label = ttk.Label(frame, width=par_name_label_width,
                               text=par_label_text)
    # Distribution menu
    if par_name == 'FigureDPI':
        pass
    else:
        distr_menu = tk.OptionMenu(
            frame, par_args_variables['distribution'],
            *distr_options, command=lambda _: controller.change_distribution(frame,
                                                                             pars_label_width=par_name_label_width))
    # Value argument label
    value_label = ttk.Label(frame, text='Value:',
                            width=distr_arg_label_width)

    if frame.par_bounds['discrete_bounds'] is not None:
        standard_tooltip_text = (
            'Set value of {}.\nPossible values are \n{}.'.format(
                frame.toolTipText,
                controller.reformat_list_presentation(
                    frame.par_bounds['discrete_bounds'])))
    else:
        standard_tooltip_text=(
            'Set value of {}.\nPossible values are between {} and {}.'.format(
                frame.toolTipText,
                frame.par_bounds['lower_bound'],
                frame.par_bounds['upper_bound']))

    # Value argument text field
    value_entry = ParameterEntry(frame, par_name, par_args_variables['value'],
                                 text_field_width, tool_tip,
                                 standard_tooltip_text=standard_tooltip_text,
                                 **frame.par_bounds)

    # Configure widgets on the frame
    if par_name == 'FigureDPI':
        pass
    else:
        distr_menu.config(width=distr_menu_width)
        tool_tip.bind(distr_menu, 'Select distribution for {}.'.format(
            frame.toolTipText))

    # Place widgets on the frame
    par_name_label.grid(row=0, column=0, padx=5, sticky='w')
    if par_name == 'FigureDPI':
        value_label.grid(row=0, column=1, padx=5, sticky='w')
        value_entry.grid(row=0, column=2, padx=5, sticky='w')
    else:
        distr_menu.grid(row=0, column=1, padx=5, sticky='w')
        value_label.grid(row=0, column=2, padx=5, sticky='w')
        value_entry.grid(row=0, column=3, padx=5, sticky='w')

    # Save frame
    controller.setvar(name='.'.join([cmpnt_nm, par_name, 'frame']), value=frame)


def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('AoR')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # Populate dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        AOR_PARAMETER_VALUES)

    for output_key in AOR_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    comp_type_label = ttk.Label(tab,
                                text="AOR Workflow Component",
                                font=LABEL_FONT,
                                name="aorTabLabel")
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(AOR_PARAMETERS):

        if par_name == 'BrineDensity' and controller.component_list['wellbore'][1].get() not in brineDensityWellbores:
            pass
        else:
            par_frames[par_name] = tk.Frame(tab, name=par_name[0].lower()+par_name[1:])
            par_frames[par_name].grid(row=ind + 1, column=0, sticky='w',
                                      padx=PARAMETER_FRAME_PADX)

            setup_parameter_frame(
                controller, par_frames[par_name], par_name,
                {'lower_bound': AOR_PARAMETER_VALUES[par_name][6],
                 'upper_bound': AOR_PARAMETER_VALUES[par_name][7]},
                AOR_PARAMETERS_SETUP[par_name][0],
                AOR_PARAMETERS_SETUP[par_name][1],
                PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
                cmpnt_nm, tool_tip)

    # Injection well control
    inj_well_frame = tk.Frame(tab, name='injWellLocation')
    inj_well_frame.grid(row=10, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Option to plot injection sites
    componentVars['plotInjectionSites'] = BooleanVar()
    componentVars['plotInjectionSites'].set(0)
    inj_site_plot_label = ttk.Label(inj_well_frame,
                                    width=PARAMETER_LABEL_WIDTH,
                                    text='Plot Injection Sites:',
                                    name='injSitePlotLabel')
    inj_site_plot_checkbox = tk.Checkbutton(inj_well_frame,
                                            variable=componentVars['plotInjectionSites'],
                                            command=lambda: disable_injection_location_widgets(controller,
                                                                                               inj_well_frame))
    inj_site_plot_label.grid(row=0, column=0, sticky='w', padx=5)
    inj_site_plot_checkbox.grid(row=0, column=1, sticky='w', padx=5)
    inj_well_frame.checkbox_variable = componentVars['plotInjectionSites']

    if controller.component_list['reservoir'][1].get() == 'Lookup Table Reservoir':
        inj_well_label_text = "Injection well location(s):"
        arg_labels = ['x-coordinate(s) [m]:', 'y-coordinate(s) [m]:']
        tool_tip_text = ''.join([
            'Enter {}-coordinate of each injection well, separated by commas.\n',
            'The number of the provided x- and y-coordinates must be the same.'])

        inj_well_label = ttk.Label(inj_well_frame, text=inj_well_label_text)
        inj_well_label.grid(row=1, column=0, sticky='w', padx=5)

        inj_well_coords_frame = ttk.Frame(inj_well_frame, name="injWellCoordsFrame")
        inj_well_coords_frame.grid(row=2, column=0, sticky='w', padx=20)

        coords = ['x', 'y']
        arg_names = ['injX', 'injY']

        coords_labels = []
        coords_fields = []
        # Create and place label and entry widgets
        for ind in range(2):
            # Create variable to keep value of the coordinate
            componentVars[cmpnt_nm][arg_names[ind]] = StringVar()
            componentVars[cmpnt_nm][arg_names[ind]].set("0")

            coords_labels.append(ttk.Label(
                inj_well_coords_frame, text=arg_labels[ind], width=LABEL_WIDTH))

            coords_fields.append(tk.Entry(
                inj_well_coords_frame,
                width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                textvariable=componentVars[cmpnt_nm][arg_names[ind]]))
            coords_labels[-1].grid(row=ind, column=0, pady=5, padx=5, sticky='w')
            coords_fields[-1].grid(row=ind, column=1, pady=5, padx=10, sticky='w')
            tool_tip.bind(coords_fields[-1],
                          tool_tip_text.format(coords[ind]))
        inj_well_frame.coords_fields = coords_fields
        disable_injection_location_widgets(controller, inj_well_frame)
