"""
Module contains several methods needed for creating tab (page) in GUI
for AOR Workflow component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import BooleanVar
from tkinter import IntVar

from openiam.gu_interface.dictionarydata import (
    componentVars, componentChoices, workflowVars, workflowChoices, connectionsDictionary, componentTypeDictionary,
    DISTRIBUTION_OPTIONS, LABEL_FONT, PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH, DISTRIBUTION_ARG_LABEL_WIDTH,
    DISTRIBUTION_ARG_TEXTFIELD_WIDTH, OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX
    )


from openiam.gu_interface.cmpnts_tabs.locations import (
    add_inj_well_frame_widgets, add_obs_locs_frame_widgets, read_obs_locations_data, read_grid_locations_data
    )
from openiam.gu_interface.cmpnts_tabs.commons import commons_read_tab_vars
from openiam.gu_interface.cmpnts_tabs.parameter_entry import ParameterEntry

LABEL_WIDTH = 17

AOR_PARAMETERS = ['FigureDPI',
                  'CriticalPressureMPa',
                  'BrineDensity',
                  'PlotInjectionSites',
                  'CriticalPressureMPa_no_calc']

AOR_PARAMETERS_SETUP = {
    'FigureDPI': ["Figure Resolution [DPI]:", "FigureDPI"],
    'CriticalPressureMPa': ["CriticalPressure [MPa]:",
                            'CriticalPressureMPa'],
    'BrineDensity': ["Brine density [kg/m{}]:".format(u'\u00B3'),
                     'Brine density'],
    }

# Set AOR Workflow parameters names and value, min, max, second value, mean, std, bounds
AOR_PARAMETER_VALUES = {
    'FigureDPI': [100, 100, 400, 100, 200, 50, 100, 400],
    'CriticalPressureMPa': [5.0, 1.0, 25.0, 10.0, 5.0, 0.01, 0.1, 90.0],
    'BrineDensity': [1000, 900, 1500, 1200, 1200, 150, 900, 1500],
    }

AOR_OBSERVATIONS = []

AOR_OBSERVATIONS_SETUP = {}

brineDensityWellbores = ['OpenWellbore', 'MultisegmentedWellbore']


def disable_injection_location_widgets(controller, frame):
    """
    Disable injection location widgets
    on the AOR Workflow Tab based on
    whether to plot injection location.
    Only applicable in the case of the
    LookupTableReservoir.
    """

    # Only control the enabling/disabling of coordinate input if
    # the Lookup Table Reservoir component is selected
    if controller.component_list['reservoir'][1].get() == 'Lookup Table Reservoir':

        # Set association between checkbutton state and frame state
        # In this case, checked means enabled and unchecked is disabled
        frame_state = {1: 'normal', 0: 'disabled'}

        # Get the current state of the checkbutton
        check_button_state = frame.checkbox_variable.get()

        # Enable or disable the widgets for input
        # of injection well x and y based on checkbutton state
        for field in frame.coords_fields:
            field.configure(state=frame_state[check_button_state])


def disable_crit_pressure_widgets(controller, frame):
    """
    Disable critical pressure entry widget
    on the AOR Workflow Tab based on
    whether "Calculate Critical Pressure"
    is checked.
    """

    # Set association between checkbutton state and frame state
    # In this case, checked means enabled and unchecked is disabled
    frame_state = {0: 'normal', 1: 'disabled'}

    # Get the current state of the checkbutton
    check_button_state = frame.checkbutton_variable.get()

    # Enable or disable the widgets for input
    # of injection well x and y based on checkbutton state
    frame.crit_press_entry.configure(state=frame_state[check_button_state])


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'AoR', parameter_names=AOR_PARAMETERS,
        observation_names=None)

    cmpnt_data = read_grid_locations_data(cmpnt_data, cmpnt_nm)

    if 'LookupTableReservoir1' in componentVars.keys():
        cmpnt_data = read_obs_locations_data(cmpnt_data, cmpnt_nm)
        cmpnt_data['InjectionLocations'] = cmpnt_data.pop('Locations')

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
    par_name_label = ttk.Label(frame, width=par_name_label_width, text=par_label_text)

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
    value_entry = ParameterEntry(frame, par_name, par_args_variables,
                                 text_field_width, tool_tip,
                                 standard_tooltip_text=standard_tooltip_text,
                                 **frame.par_bounds)


    # Place widgets on the frame
    par_name_label.grid(row=0, column=0, padx=5, sticky='w')
    value_entry.grid(row=0, column=1, padx=5, sticky='w')

    # Save frame
    controller.setvar(name='.'.join([cmpnt_nm, par_name, 'frame']), value=frame)


def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append(cmpnt_type)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    #Populate dictionary
    componentVars[cmpnt_nm]['Params'] = {}

    comp_type_label = ttk.Label(tab,
                                text="AOR Workflow Components",
                                font=LABEL_FONT,
                                name="aorTabLabel")
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}

    # Count row of current frame
    ind = 1

    # FigureDPI Entry
    componentVars[cmpnt_nm]['Params']['FigureDPI'] = IntVar()
    componentVars[cmpnt_nm]['Params']['FigureDPI'].set(AOR_PARAMETER_VALUES['FigureDPI'][0])
    par_frames['FigureDPI'] = tk.Frame(tab, name='figureDPI')
    par_frames['FigureDPI'].grid(row=ind, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    setup_parameter_frame(controller, par_frames['FigureDPI'], 'FigureDPI',
                {'lower_bound': AOR_PARAMETER_VALUES['FigureDPI'][6],
                          'upper_bound': AOR_PARAMETER_VALUES['FigureDPI'][7]},
                          AOR_PARAMETERS_SETUP['FigureDPI'][0],
                          AOR_PARAMETERS_SETUP['FigureDPI'][1],
                          PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                          DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                          DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params']['FigureDPI'],
                          cmpnt_nm, tool_tip)
    ind += 1

    # Critical Pressure Entry
    componentVars[cmpnt_nm]['Params']['CriticalPressureMPa'] = BooleanVar()
    componentVars[cmpnt_nm]['Params']['CriticalPressureMPa'].set(1)
    componentVars[cmpnt_nm]['Params']['CriticalPressureMPa_no_calc'] = StringVar()
    componentVars[cmpnt_nm]['Params']['CriticalPressureMPa_no_calc'].set(AOR_PARAMETER_VALUES['CriticalPressureMPa'][0])
    crit_press_tooltip_text = ('Set value of {}.\nPossible values are between {} and {}.'.format(
                         AOR_PARAMETERS_SETUP['CriticalPressureMPa'][1],
                               AOR_PARAMETER_VALUES['CriticalPressureMPa'][1],
                               AOR_PARAMETER_VALUES['CriticalPressureMPa'][2]))
    crit_press_bounds = {'lower_bound': AOR_PARAMETER_VALUES['CriticalPressureMPa'][1],
                         'upper_bound': AOR_PARAMETER_VALUES['CriticalPressureMPa'][2],
                         'discrete_bounds': None}

    par_frames['CriticalPressureMPa'] = tk.Frame(tab, name='criticalPressureMPa')
    par_frames['CriticalPressureMPa'].grid(row=ind, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    crit_press_label1 = ttk.Label(par_frames['CriticalPressureMPa'],
                                 width=PARAMETER_LABEL_WIDTH,
                                 text=AOR_PARAMETERS_SETUP['CriticalPressureMPa'][0])

    crit_press_checkbox = tk.Checkbutton(par_frames['CriticalPressureMPa'],
                                         variable=componentVars[cmpnt_nm]['Params']['CriticalPressureMPa'],
                                         command=lambda: disable_crit_pressure_widgets(controller,
                                                                                       par_frames['CriticalPressureMPa']),
                                         name="crit_press_check")

    crit_press_label2 = ttk.Label(par_frames['CriticalPressureMPa'],
                                 width=PARAMETER_LABEL_WIDTH,
                                 text='Calculated')

    crit_press_entry = ParameterEntry(par_frames['CriticalPressureMPa'], 'CriticalPressureMPa',
                                      componentVars[cmpnt_nm]['Params']['CriticalPressureMPa_no_calc'],
                                      DISTRIBUTION_ARG_TEXTFIELD_WIDTH, tool_tip,
                                      standard_tooltip_text=crit_press_tooltip_text,
                                      **crit_press_bounds)

    par_frames['CriticalPressureMPa'].checkbutton_variable = componentVars[cmpnt_nm]['Params']['CriticalPressureMPa']
    par_frames['CriticalPressureMPa'].crit_press_entry = crit_press_entry
    disable_crit_pressure_widgets(controller, par_frames['CriticalPressureMPa'])

    # Place widgets on the frame
    crit_press_label1.grid(row=0, column=0, padx=5, sticky='w')
    crit_press_checkbox.grid(row=0, column=1, padx=5, sticky='w')
    crit_press_label2.grid(row=0, column=2, padx=5, sticky='w')
    crit_press_entry.grid(row=0, column=3, padx=5, sticky='w')

    # Save frame
    controller.setvar(name='.'.join([cmpnt_nm, 'CriticalPressureMPa', 'frame']), value=par_frames['CriticalPressureMPa'])

    # Brine Density Entry
    if controller.component_list['reservoir'][1].get().replace(" ", "") == 'AnalyticalReservoir' and \
        controller.component_list['wellbore'][1].get().replace(" ", "") in brineDensityWellbores:

        ind += 1
        componentVars[cmpnt_nm]['Params']['BrineDensity'] = StringVar()
        componentVars[cmpnt_nm]['Params']['BrineDensity'].set(AOR_PARAMETER_VALUES['BrineDensity'][0])
        brine_density_tooltip_text = ('Set value of {}.\nPossible values are between {} and {}.'.format(
                                AOR_PARAMETERS_SETUP['BrineDensity'][1],
                                      AOR_PARAMETER_VALUES['BrineDensity'][1],
                                      AOR_PARAMETER_VALUES['BrineDensity'][2]))
        brine_density_bounds = {'lower_bound': AOR_PARAMETER_VALUES['BrineDensity'][1],
                                'upper_bound': AOR_PARAMETER_VALUES['BrineDensity'][2],
                                'discrete_bounds': None}

        par_frames['BrineDensity'] = tk.Frame(tab, name='brineDensity')
        par_frames['BrineDensity'].grid(row=ind, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

        brine_density_label = ttk.Label(par_frames['BrineDensity'],
                                        width=PARAMETER_LABEL_WIDTH,
                                        text=AOR_PARAMETERS_SETUP['BrineDensity'][0])

        brine_density_entry = ParameterEntry(par_frames['BrineDensity'], 'BrineDensity',
                                             componentVars[cmpnt_nm]['Params']['BrineDensity'],
                                             DISTRIBUTION_ARG_TEXTFIELD_WIDTH, tool_tip,
                                             standard_tooltip_text=brine_density_tooltip_text,
                                             **brine_density_bounds)

        # Place widgets on the frame
        brine_density_label.grid(row=0, column=0, padx=5, sticky='w')
        brine_density_entry.grid(row=0, column=1, padx=5, sticky='w')

        # Save frame
        controller.setvar(name='.'.join([cmpnt_nm, 'BrineDensity', 'frame']),
                          value=par_frames['BrineDensity'])

    # Wellbore Location
    ind += 1
    componentVars[cmpnt_nm]['Locations'] = {'grid': {}}
    well_loc_frame = tk.Frame(tab, name='wellboreLocation')
    well_loc_frame.grid(row=ind, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    well_loc_label_text = "Wellbore grid locations:"
    tool_tip_text = ''.join([
        'Enter the {} of the wellbore grid in {}.',
        'WARNING: If using the Generic Aquifer component, the simulation will run very',
        'slowly if grid spacing is small (large number of wellbores simulated).'])

    well_loc_label = ttk.Label(well_loc_frame, text=well_loc_label_text)
    well_loc_label.grid(row=0, column=0, sticky='w', padx=5)

    well_loc_coords_frame = ttk.Frame(well_loc_frame, name="wellboreCoordsFrame")
    well_loc_coords_frame.grid(row=1, column=0, sticky='w', padx=20)

    coords = ['x', 'y']
    loc_type = ['minimum', 'maximum', 'number of wells']
    arg_labels = ['x-min [m]:', 'x-max [m]:', 'x size:', 'y-min [m]:', 'y-max [m]:', 'y size:']
    arg_names = ['xmin', 'xmax', 'xsize', 'ymin', 'ymax', 'ysize']

    coords_labels = []
    coords_fields = []
    # Create and place label and entry widgets
    for r in range(2):
        for l in range(3):
            if r == 0:
                n = l
            else:
                n = l + 3
            # Create variable to keep value of the coordinate
            if arg_names[n][-4:] != 'size':
                componentVars[cmpnt_nm]['Locations']['grid'][arg_names[n]] = StringVar()
                if arg_names[n][-3:] == 'min':
                    componentVars[cmpnt_nm]['Locations']['grid'][arg_names[n]].set("-50")
                else:
                    componentVars[cmpnt_nm]['Locations']['grid'][arg_names[n]].set("50")
            else:
                componentVars[cmpnt_nm]['Locations']['grid'][arg_names[n]] = IntVar()
                componentVars[cmpnt_nm]['Locations']['grid'][arg_names[n]].set(2)

            coords_labels.append(ttk.Label(
                well_loc_coords_frame, text=arg_labels[n], width=LABEL_WIDTH))

            coords_fields.append(tk.Entry(
                well_loc_coords_frame,
                width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                textvariable=componentVars[cmpnt_nm]['Locations']['grid'][arg_names[n]]))
            coords_labels[-1].grid(row=n, column=0, pady=5, padx=40, sticky='e')
            coords_fields[-1].grid(row=n, column=1, pady=5, sticky='w')
            tool_tip.bind(coords_fields[-1],
                          tool_tip_text.format(loc_type[l], coords[r]))
    well_loc_frame.coords_fields = coords_fields

    # Injection well control
    ind += 1
    inj_well_frame = tk.Frame(tab, name='injWellLocation')
    inj_well_frame.grid(row=ind, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Option to plot injection sites
    componentVars[cmpnt_nm]['Params']['PlotInjectionSites'] = BooleanVar()
    componentVars[cmpnt_nm]['Params']['PlotInjectionSites'].set(0)
    inj_site_plot_label = ttk.Label(inj_well_frame,
                                    width=PARAMETER_LABEL_WIDTH,
                                    text='Plot Injection Sites:',
                                    name='injSitePlotLabel')
    inj_site_plot_checkbox = tk.Checkbutton(inj_well_frame,
                                            variable=componentVars[cmpnt_nm]['Params']['PlotInjectionSites'],
                                            command=lambda: disable_injection_location_widgets(controller,
                                                                                               inj_well_frame))
    inj_site_plot_label.grid(row=0, column=0, sticky='w', padx=5)
    inj_site_plot_checkbox.grid(row=0, column=1, sticky='w', padx=5)
    inj_well_frame.checkbox_variable = componentVars[cmpnt_nm]['Params']['PlotInjectionSites']

    if controller.component_list['reservoir'][1].get().replace(" ", "") == 'LookupTableReservoir':
        inj_well_label_text = "Injection well location(s):"
        arg_labels = ['x-coordinate(s) [m]:', 'y-coordinate(s) [m]:']
        tool_tip_text = ''.join([
            'Enter {}-coordinate of each injection well to plot, separated by commas.\n',
            'The number of the provided x- and y-coordinates must be the same.'])

        inj_well_label = ttk.Label(inj_well_frame, text=inj_well_label_text)
        inj_well_label.grid(row=1, column=0, sticky='w', padx=5)

        inj_well_coords_frame = ttk.Frame(inj_well_frame, name="injWellCoordsFrame")
        inj_well_coords_frame.grid(row=2, column=0, sticky='w', padx=20)

        coords = ['x', 'y']
        arg_names = ['xCoordinates', 'yCoordinates']

        coords_labels = []
        coords_fields = []
        # Create and place label and entry widgets
        for r in range(2):
            # Create variable to keep value of the coordinate
            componentVars[cmpnt_nm][arg_names[r]] = StringVar()
            componentVars[cmpnt_nm][arg_names[r]].set("0")

            coords_labels.append(ttk.Label(
                inj_well_coords_frame, text=arg_labels[r], width=LABEL_WIDTH))

            coords_fields.append(tk.Entry(
                inj_well_coords_frame,
                width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                textvariable=componentVars[cmpnt_nm][arg_names[r]]))
            coords_labels[-1].grid(row=r, column=0, pady=5, padx=40, sticky='e')
            coords_fields[-1].grid(row=r, column=1, pady=5, sticky='w')
            tool_tip.bind(coords_fields[-1],
                          tool_tip_text.format(coords[r]))
        inj_well_frame.coords_fields = coords_fields
        disable_injection_location_widgets(controller, inj_well_frame)
