"""
THIS SCRIPT SERVES AS A TEMPLATE FOR CREATING WORKFLOW TABS TO BE USED
IN OPENIAM WORKFLOWS. THIS TEMPLATE PROVIDES THE BASIC REQUIRED FEATURES
BUT MAY NOT COVER ALL CASES. ADDITIONAL CODING MAY BE NECESSARY TO FIT
YOUR PARTICULAR WORKFLOW USE CASE.
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

# SET THE NAMES OF THE PARAMETERS TO BE USED IN THE WORKFLOW TAB
# THESE PARAMETER NAMES WILL BE USED AS KEYS IN THE DICTIONARY
# THAT IS OUTPUT TO THE YAML AND OPENIAM FILES IN THE WORKFLOW
# OPTIONS SECTION, SO ENSURE THAT THEY MATCH THE EXPECTED INPUTS
# FOR THE WORKFLOW.PY FILE
# REPLACE "WORKFLOW" WITH THE THREE LETTER WORKFLOW ACRONYM
# REPLACE 'parameter1' AND 'parameter2' WITH THE PARAMETERS AND
# ADD ADDITIONAL PARAMETERS AS NECESSARY
WORKFLOW_PARAMETERS = ['parameter1',
                       'parameter2']

# SET THE LABEL AND TOOL TIP NAME FOR THE PARAMETERS IN WORKFLOW_PARAMETERS
# REPLACE 'parameter1' AND 'parameter2' WITH THE PARAMETERS,
# REPLACE 'Parameter 1 Name' AND 'Parameter 2 Name' WITH THE LABELS FOR PARAMETER INPUTS,
# AND REPLACE "Parameter 1" AND "Parameter 2" WITH THE DISPLAY NAME OF THE PARAMETERS
# FOR THE TOOL TIP

WORKFLOW_PARAMETERS_SETUP = {
    'parameter1': ["Parameter 1 Name",
                   "Parameter 1"],
    'parameter2': ["Parameter 2 Name",
                   "Parameter 2"],
    }

# SET THE DEFAULT VALUE AND ACCEPTABLE BOUNDS OF EACH PARAMETER
# WORKFLOWS EXPECT A SINGLE VALUE WITH MINIMUM AND MAXIMUM VALUES
# UNLESS A SPECIALIZED EXCEPTION IS MADE
# REPLACE "WORKFLOW" WITH THE THREE LETTER WORKFLOW ACRONYM
# REPLACE 'parameter1' AND 'parameter2' WITH THE PARAMETERS
# FOR EACH PARAMETER, REPLACE 'p#_value', 'p#_min' AND 'p#_max' WITH
# THE CORRESPONDING DEFAULT VALUE, MINIMUM AND MAXIMUM VALUES
WORKFLOW_PARAMETER_VALUES = {
    'parameter1': [p1_value, p1_min, p1_max],
    'parameter2': [p2_value, p2_min, p2_max],
    }

# DEFINE ANY OTHER VARIABLES OR SETS OF VARIABLES NECESSARY HERE
# THAT NEED TO BE ACCESSED BY ALL FUNCTIONS

# THIS IS AN EXAMPLE OF TKINTER CONTROL CODE FOR ENABLING/DISABLING
# WIDGETS BASED ON THE TYPE OF COMPONENT(S) SELECTED FOR THE WORKFLOW.
# WHILE THIS FUNCTION MAY NOT BE NECESSARY FOR EVERY WORKFLOW, IT IS
# USEFUL TO NOTE THAT THIS FUNCTIONALITY IS AVAILABLE.
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


def read_tab_vars(cmpnt_nm):

    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'AoR', parameter_names=WORKFLOW_PARAMETERS,
        observation_names=None)

    return cmpnt_data


# THIS FUNCTION ALLOWS FOR THE LINE-BY-LINE SETUP OF PARAMETERS
# TO GATHER USER INPUT FOR THE WORKFLOW PARAMETERS
def setup_parameter_frame(controller, frame, par_name, par_bounds,
                          par_label_text, tool_tip_text,
                          par_name_label_width, distr_menu_width,
                          distr_arg_label_width, text_field_width,
                          distr_options, par_args_variables, cmpnt_nm, tool_tip):
    """ Create, configure and place widgets for a given parameter frame."""
    # Update attributes of frame

    # Label text for parameter
    frame.labelText = par_label_text

    # Corresponding component for parameter
    frame.component = cmpnt_nm

    # Distribution type for parameter (if necessary)
    frame.distType = par_args_variables['distribution']

    # Parameters of current component
    frame.paramVars = par_args_variables

    # Tooltip text that defines the parameter being set
    frame.toolTipText = tool_tip_text

    # Allowable upper, lower, and discrete bounds of parameter input
    frame.par_bounds = {'lower_bound': None,
                        'upper_bound': None,
                        'discrete_bounds': None}

    # Available distribution options for parameter (if necessary)
    frame.distr_options = distr_options

    # Copy provided values for upper, lower, and discrete bounds
    for key in par_bounds:
        frame.par_bounds[key] = par_bounds[key]

    # Create common widgets
    # Label with parameter name and units
    par_name_label = ttk.Label(frame, width=par_name_label_width,
                               text=par_label_text)

    # If discrete bounds is used, use this tooltip text
    if frame.par_bounds['discrete_bounds'] is not None:
        standard_tooltip_text = (
            'Set value of {}.\nPossible values are \n{}.'.format(
                frame.toolTipText,
                controller.reformat_list_presentation(
                    frame.par_bounds['discrete_bounds'])))

    # Otherwise, use the tooltip text with upper and lower bounds
    else:
        standard_tooltip_text=(
            'Set value of {}.\nPossible values are between {} and {}.'.format(
                frame.toolTipText,
                frame.par_bounds['lower_bound'],
                frame.par_bounds['upper_bound']))

    # Sets up entry field for parameter value including validation using upper, lower and discrete bounds
    value_entry = ParameterEntry(frame, par_name, par_args_variables['value'],
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

    # tab: address of workflow tab within the larger Tkinter tab/window set
    # cmpnt_nm: Workflow
    # cmpnt_type: three-letter acronym for the current workflow type
    # tool_tip: Tooltip text that defines the parameter being set
    # args: Additional parameters that may be necessary for the workflow setup

    # Add the Workflow component to the componentChoices list
    componentChoices.append(cmpnt_nm)

    # Add the workflow type
    componentTypeDictionary.append(cmpnt_type)

    # Add a key to the componentVars dictionary to store the parameters
    # associated with the workflow for later output
    componentVars[cmpnt_nm] = {}

    # Add the component name and component type to the componentVars dictionary
    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # Add the parameters to be set in the Workflow tab to the componentVars dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
       WORKFLOW_PARAMETER_VALUES)

    # Create and place the main label for the workflow tab
    comp_type_label = ttk.Label(tab,
                                text=cmpnt_nm + " Workflow Component",
                                font=LABEL_FONT,
                                name=cmpnt_nm.lower() + "TabLabel")
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    # Create dictionary to hold the frames for each parameter of the workflow
    par_frames = {}

    # Create a frame for each parameter of the workflow and add it to the dictionary
    for ind, par_name in enumerate(WORKFLOW_PARAMETERS):

        par_frames[par_name] = tk.Frame(tab, name=par_name[0].lower()+par_name[1:])
        par_frames[par_name].grid(row=ind + 1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        setup_parameter_frame(
            controller, par_frames[par_name], par_name,
            {'lower_bound': WORKFLOW_PARAMETER_VALUES[par_name][6],
             'upper_bound': WORKFLOW_PARAMETER_VALUES[par_name][7]},
            WORKFLOW_PARAMETERS_SETUP[par_name][0],
            WORKFLOW_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Any additional parameter settings not covered in the
    # setup_parameter_frame function can be handled here
    # Included below is an example of controlling whether
    # to plot injection sites and turning on/off the input
    # for well locations based on the checkbox state

    # Injection well control
    # Create and place label for the injection well location frame
    inj_well_frame = tk.Frame(tab, name='injWellLocation')
    inj_well_frame.grid(row=10, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Create and set the Boolean variable control the checkbox state
    componentVars['plotInjectionSites'] = BooleanVar()
    componentVars['plotInjectionSites'].set(0)

    # Create and place the label and checkbox for plotting injection site
    # The "disable_injection_widgets" function used in the checkbox "command"
    # value turns on/off the subsequent entries for well location based
    # on the current state of the checkbox Boolean variable
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

    # This is an example of component-specific option creation
    # Only create x and y coordinate inputs if the reservoir chosen is Lookup Table
    if controller.component_list['reservoir'][1].get() == 'Lookup Table Reservoir':

        # Create and place the label for injection well location input
        inj_well_label_text = "Injection well location(s):"
        arg_labels = ['x-coordinate(s) [m]:', 'y-coordinate(s) [m]:']
        tool_tip_text = ''.join([
            'Enter {}-coordinate of each injection well, separated by commas.\n',
            'The number of the provided x- and y-coordinates must be the same.'])

        inj_well_label = ttk.Label(inj_well_frame, text=inj_well_label_text)
        inj_well_label.grid(row=1, column=0, sticky='w', padx=5)

        # Create a sub-frame to hold the x and y coordinate input
        inj_well_coords_frame = ttk.Frame(inj_well_frame, name="injWellCoordsFrame")
        inj_well_coords_frame.grid(row=2, column=0, sticky='w', padx=20)

        # Create lists of names for x and y coordinate labels and variables
        coords = ['x', 'y']
        arg_names = ['injX', 'injY']

        # Create lists to hold the labels and inputs for x and y coordinates
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

        # Add the coordinate fields to the parent frame
        inj_well_frame.coords_fields = coords_fields

        # Run the code to disable these inputs to correspond with the
        # "Plot Injection Sites" checkbox starting as unchecked
        disable_injection_location_widgets(controller, inj_well_frame)
