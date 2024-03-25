# Pmw is imported to enable hover text
# pickle is used to generate and read binary files
import os
import sys
import pickle
import random
import logging
import warnings
import yaml
logging.basicConfig(level=logging.WARNING)

import tkinter as tk
from tkinter import ttk, StringVar, IntVar, DoubleVar, BooleanVar, messagebox

import matplotlib
debug_msg = 'Available matplotlib backends: {}'.format(matplotlib.rcsetup.all_backends)
logging.debug(debug_msg)

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    # matplotlib.use('TkAgg')

import numpy as np
import Pmw

from openiam.gu_interface.Disclaimer import Disclaimer_Page
from openiam.gu_interface.Dashboard import Dashboard_Page
from openiam.gu_interface.OpenIAM_Page import OpenIAM_Page, disable_time_frame_widgets
from openiam.gu_interface.PostProcessor_Page import PostProcessor_Page
from openiam.gu_interface.Workflow_Page import Workflow_Page

from openiam.gu_interface.dictionarydata import (
    d, APP_SIZE, TAB_SIZE, componentVars, componentChoices,
    componentTypeDictionary, connectionsDictionary,
    DISTRIBUTION_OPTIONS, connections, connectionTypes,
    COMPONENT_TYPES, ANALYSIS_TYPES,
    DISTRIBUTION_MENU_WIDTH, DISTRIBUTION_ARG_LABEL_WIDTH,
    DISTRIBUTION_ARG_TEXTFIELD_WIDTH, PARAMETER_LABEL_WIDTH,
    GFR_PARAMETER_LABEL_WIDTH, STRATA_PARAMETER_LABEL_WIDTH,
    FL_PARAMETER_LABEL_WIDTH,
    MODEL_TAB_LABEL_WIDTH2, MODEL_TAB_LABEL_WIDTH3,
    MODEL_TAB_ENTRY_WIDTH, MODEL_TAB_MENU_WIDTH,
    DISTRIBUTION_PARS_LABELS, DISTRIBUTION_PARS_SETUPS,
    workflowVars, WORKFLOW_TYPES)

from openiam.gu_interface.cmpnts_tabs import (
    src_tab, arc_tab, grc_tab, trc_tab, msw_tab, lutr_tab, cw_tab,
    cwwr_tab, ow_tab, gfr_tab, ff_tab, fl_tab, hcl_tab, sh_tab, ca_tab,
    aalf_tab, daa_tab, daaml_tab, fgaq_tab, fgaz_tab, ga_tab, atm_tab, 
    psa_tab, strata_tab, cws_tab, AOR_wf_tab, locations, wf_strata_tab)
from openiam.gu_interface.cmpnts_tabs.parameter_entry import ParameterEntry

from openiam.components.iam_base_classes import IAM_DIR

# Save location of source folder in the top level folder
USER_DIR = os.sep.join([IAM_DIR, 'examples', 'user'])

# Dictionaries to set aquifer output based on component and plume type
TTFD_AQUIFER_COMPONENT_OUTPUT = {
    'FutureGen2Aquifer': {'pH': ['pH_dx', 'pH_dy', 'pH_dz'],
                          'TDS': ['TDS_dx', 'TDS_dy', 'TDS_dz'],
                          'Dissolved_CO2': ['Dissolved_CO2_dx',
                                            'Dissolved_CO2_dy',
                                            'Dissolved_CO2_dz'],
                          'Pressure': ['Pressure_dx', 'Pressure_dy', 'Pressure_dz']},
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

def set_aquifer_outputs(value, controller):
    """
    Write a function to set the aquifer
    component outputs based on the
    aquifer type and the plume type.
    """
    aquifer_type = controller.component_list['aquifer'][1].get().replace(" ", "")
    aq_plume_selection = value
    aq_output_types = TTFD_AQUIFER_COMPONENT_OUTPUT[aquifer_type]
    aq_output_types = aq_output_types[aq_plume_selection]

    # Reset outputs of aquifer based on aquifer type
    if aquifer_type in ["CarbonateAquifer", "DeepAlluviumAquifer"]:
        aq_output_num = 12
    elif aquifer_type == 'GenericAquifer':
        aq_output_num = 8
    elif aquifer_type == 'FutureGen2Aquifer':
        aq_output_num = 16
    elif aquifer_type == 'FutureGen2AZMI':
        aq_output_num = 20

    aquifer_checkbuttons = []
    aquifer_checkframe = controller.workflow_tabs[2].winfo_children()[0] \
                          .winfo_children()[0].winfo_children()[-1].winfo_children()

    for f in aquifer_checkframe:
        if type(f) == tk.Checkbutton:
            aquifer_checkbuttons.append(f)

    aq_output_keys = list(componentVars[aquifer_type + "1"].keys())[-aq_output_num+1:-1]
    aq_dict = {k: v for k, v in zip(aq_output_keys, aquifer_checkbuttons)}

    # Select new outputs
    for curr_key in aq_output_keys:
        if curr_key in aq_output_types:
            aq_dict[curr_key].config(state='normal')
            componentVars[aquifer_type+"1"][curr_key].set(1)
            aq_dict[curr_key].config(state='disabled')
        else:
            aq_dict[curr_key].config(state='normal')
            componentVars[aquifer_type+"1"][curr_key].set(0)
            aq_dict[curr_key].config(state='disabled')

def workflow_init(controller, data):
    component_list = {}
    aquiferName = StringVar()
    ckeys = data['ModelParams']['Components'] + ['Workflow']
    for key in ckeys:
        if key == "Workflow":
            data[key]['type'] = data[key].pop('Type')

        if aquiferName.get() == '' or aquiferName.get() == 'none':
            try:
                aquiferName.set(data[key]['AquiferName'])
            except KeyError:
                try:
                    aquiferName.set(data[key]['LeakTo'])
                except KeyError:
                    aquiferName.set('none')

        if key.find('Reservoir') != -1:
            comp_name = 'reservoir'
        elif key.find('Wellbore') != -1:
            comp_name = 'wellbore'
        elif key.find('Aquifer') != -1 or key.find('AZMI') != -1:
            comp_name = 'aquifer'
        elif key.find('Workflow') != -1:
            comp_name = 'Workflow'

        component_list[comp_name] = [StringVar(), StringVar()]
        component_list[comp_name][0].set(key)
        component_list[comp_name][1].set(data[key]['type'])

    conn = StringVar()
    conn.set('Auto')

    controller.add_workflow(conn, aquiferName, controller.wftabControl,
                            controller.wftabControl.connection_menu,
                            controller.wftabControl.workflowSetupFrame, controller,
                            [], {}, component_list)

    for key in ckeys:
        # Call additional widgets setup for selected components
        if data[key]['type'] in ['AnalyticalReservoir', 'LookupTableReservoir']:
            locations.load_obs_locations_data(data[key], key)

            if data[key]['type'] == 'LookupTableReservoir':
                code, msg = lutr_tab.finish_load_setup(controller, data[key], key)
                if code == 0:
                    controller.frames[Workflow_Page].tkraise()
                    messagebox.showerror("Error", msg)
                    break
                continue

        # For compatibility with older versions of GUI files before change of the
        # parameter name brineResSaturation to aquBrineResSaturation
        if data[key]['type'] == 'MultisegmentedWellbore':
            if 'brineResSaturation' in data[key]['Parameters']:
                if isinstance(data[key]['Parameters']['brineResSaturation'], dict):
                    data[key]['Parameters']['aquBrineResSaturation'] = \
                        data[key]['Parameters']['brineResSaturation'].copy()
                else:
                    data[key]['Parameters']['aquBrineResSaturation'] = \
                        data[key]['Parameters']['brineResSaturation']
                data[key]['Parameters'].pop('brineResSaturation')


def load_workflow(controller, data):
    ckeys = data['ModelParams']['Components'] + ['Workflow']
    resvar = controller.component_list['reservoir'][0].get()
    wellvar = controller.component_list['wellbore'][0].get()

    if data['Workflow']['type'] == 'AoR' or data['Workflow']['type'] == 'TTFD':

        for key in ckeys:
            # Set up Workflow tab

            # Load figure resolution
            if key == 'Workflow':
                # Set figure DPI
                componentVars[key]['Params']['FigureDPI'].set(data[key]['Options']['FigureDPI'])

                # Set critical pressure strategy
                if data[key]['Options']['CriticalPressureMPa'] == 'Calculated':
                    componentVars[key]['Params']['CriticalPressureMPa'].set(True)
                    if 'OpenWellbore1' in ckeys:
                        get_widgets = controller.wftabControl.nametowidget(''.join(['.!frame.!workflow_page.',
                                                                                    'workflow_notebook.',
                                                                                    'openwellbore1.',
                                                                                    'openwellbore1_canvas.',
                                                                                    'openwellbore1_tabType']))
                        get_widgets.winfo_children()[8].winfo_children()[1].config(state="normal")
                        componentVars['OpenWellbore1']['Controls']['critPressureApproach'].set(1)
                        get_widgets.winfo_children()[8].winfo_children()[1].config(state="disabled")


                else:
                    componentVars[key]['Params']['CriticalPressureMPa'].set(False)
                    crit_press = float(data[key]['Options']['CriticalPressureMPa'])
                    cp_widgets = controller.wftabControl.nametowidget(''.join(['.!frame.!workflow_page.',
                                                                                'workflow_notebook.',
                                                                                'workflow.',
                                                                                'workflow_canvas.',
                                                                                'workflow_tabType.',
                                                                                'criticalPressureMPa']))
                    cp_widgets.winfo_children()[-1].config(state="normal")
                    componentVars[key]['Params']['CriticalPressureMPa_no_calc'].set(crit_press)
                    if 'OpenWellbore1' in ckeys:
                        get_widgets = controller.wftabControl.nametowidget(''.join(['.!frame.!workflow_page.',
                                                                                    'workflow_notebook.',
                                                                                    'openwellbore1.',
                                                                                    'openwellbore1_canvas.',
                                                                                    'openwellbore1_tabType']))

                        get_widgets.winfo_children()[7].winfo_children()[-1].config(state="normal")
                        componentVars['OpenWellbore1']['Params']['critPressure']['value'].set(crit_press * (10 ** 6))
                        get_widgets.winfo_children()[7].winfo_children()[-1].config(state="disabled")

                        get_widgets.winfo_children()[8].winfo_children()[1].config(state="normal")
                        componentVars['OpenWellbore1']['Controls']['critPressureApproach'].set(1)
                        get_widgets.winfo_children()[8].winfo_children()[1].config(state="disabled")

                        get_widgets.winfo_children()[9].winfo_children()[1].config(state="normal")
                        componentVars['OpenWellbore1']['Controls']['enforceCritPressure'].set(1)
                        get_widgets.winfo_children()[9].winfo_children()[1].config(state="disabled")

                if 'AnalyticalReservoir1' in ckeys:
                    brine_density = data['AnalyticalReservoir1']['Parameters']['brineDensity']['value']
                    componentVars['AnalyticalReservoir1']['Params']['brineDensity']['value'].set(brine_density)
                    if 'OpenWellbore1' in ckeys:
                        get_widgets = controller.wftabControl.nametowidget(''.join(['.!frame.!workflow_page.',
                                                                                    'workflow_notebook.',
                                                                                    'openwellbore1.',
                                                                                    'openwellbore1_canvas.',
                                                                                    'openwellbore1_tabType']))

                        get_widgets.winfo_children()[5].winfo_children()[-1].config(state="normal")
                        componentVars['OpenWellbore1']['Params']['brineDensity']['value'].set(brine_density)


                    elif 'MultisegmentedWellbore1' in ckeys:
                        get_widgets = controller.wftabControl.nametowidget(''.join(['.!frame.!workflow_page.',
                                                                                    'workflow_notebook.',
                                                                                    'multisegmentedwellbore1.',
                                                                                    'multisegmentedwellbore1_canvas.',
                                                                                    'multisegmentedwellbore1_tabType']))

                        get_widgets.winfo_children()[4].winfo_children()[-1].config(state="normal")
                        componentVars['MultisegmentedWellbore1']['Params']['brineDensity']['value'].set(brine_density)

                if data['Workflow']['type'] == 'TTFD':
                    # set plume type
                    componentVars[key]['Params']['PlumeType'].set(data[key]['Options']['PlumeType'])
                    set_aquifer_outputs(data[key]['Options']['PlumeType'], controller)

                    # set monitoring locations
                    for name in ['coordx', 'coordy', 'coordz', 'HorizontalWindow', 'VerticalWindow']:
                        temp_monitoring_loc = str(data[key]['Options']['MonitoringLocations'][name])
                        componentVars[key]['Params']['MonitoringLocations'][name].set(temp_monitoring_loc)

                    # set number of z-points in aquifers
                    componentVars[key]['Params']['NumZPointsWithinAquifers'].set(
                        data[key]['Options']['NumZPointsWithinAquifers']
                    )

                    # set number of z-points in shales
                    componentVars[key]['Params']['NumZPointsWithinShales'].set(
                        data[key]['Options']['NumZPointsWithinShales']
                    )

                    # set option for saving CSV files
                    componentVars[key]['Params']['SaveCSVFiles'].set(
                        data[key]['Options']['SaveCSVFiles']
                    )

                    # set option for writing Dream output
                    componentVars[key]['Params']['WriteDreamOutput'].set(
                        data[key]['Options']['WriteDreamOutput']
                    )

            # Load injection location information
            grid_names = ['xmin', 'xmax', 'xsize', 'ymin', 'ymax', 'ysize']
            for n in grid_names:
                componentVars['Workflow']['Locations']['grid'][n].set(data[wellvar]['Locations']['grid'][n])

            get_widgets = controller.wftabControl.nametowidget(''.join(['.!frame.!workflow_page.',
                                                                        'workflow_notebook.',
                                                                        wellvar.lower(), '.',
                                                                        wellvar.lower(), '_canvas.',
                                                                        wellvar.lower(), '_tabType']))

            if wellvar in ['openwellbore1', 'multisegmentedwellbore1']:
                wellvar_num = 12
            else:
                wellvar_num = 6

            get_widgets.winfo_children()[wellvar_num].winfo_children()[-1].config(state="normal")
            componentVars[wellvar]['number'].set(int(data[wellvar]['Locations']['grid']['xsize'] *
                                                 data[wellvar]['Locations']['grid']['ysize']))
            get_widgets.winfo_children()[wellvar_num].winfo_children()[-1].config(state="disabled")

            # Set value for plotting injection sites
            plotInjSites = data['Workflow']['Options']['PlotInjectionSites']
            componentVars['Workflow']['Params']['PlotInjectionSites'].set(plotInjSites)

            if resvar == 'LookupTableReservoir1':
                inj_names = ['InjectionCoordx', 'InjectionCoordy']
                coord_names = ['xCoordinates', 'yCoordinates']
                for n in range(len(coord_names)):
                    componentVars['Workflow'][coord_names[n]].set(data['Workflow']['Options'][inj_names[n]])

    if data['Workflow']['type'] == 'AoR':
        controller.workflowType.set('Area of Review')
    elif data['Workflow']['type'] == 'TTFD':
        controller.workflowType.set('Time to First Detection')

    add_wf_button = controller.wftabControl.nametowidget(''.join(['.!frame.!workflow_page.',
                                                                  'workflow_notebook.',
                                                                  'workflow_tab.',
                                                                  'addWorkflow_frame.',
                                                                  'addWorkflow_button']))
    add_wf_button.config(state="disabled")

    try:
        componentVars['wf_timePointsInput'].set(data['ModelParams']['TimePointsInput'])
        componentVars['wf_timePoints'].set(data['ModelParams']['TimePoints'])
    except:
        componentVars['wf_endTime'].set(data['ModelParams']['EndTime'])
        componentVars['wf_timeStep'].set(data['ModelParams']['TimeStep'])
        componentVars['wf_timePointsInput'].set(False)
        componentVars['wf_timePoints'].set('')

    componentVars['wf_analysis']['wf_type'].set(data['ModelParams']['Analysis']['type'])
    componentVars['wf_logging'].set(data['ModelParams']['Logging'])
    componentVars['wf_outputDirectory'].set(data['ModelParams']['OutputDirectory'])
    try:
        componentVars['wf_outputDirectoryGenerate'].set(data['ModelParams']['OutputDirectoryGenerate'])
    except:
        componentVars['wf_outputDirectoryGenerate'].set(False)
    controller.wf_OutputType.set(data['ModelParams']['OutputType'])
    controller.wf_GenerateOutputFiles.set(data['ModelParams']['GenerateOutputFiles'])
    controller.wf_GenerateCombOutputFile.set(data['ModelParams']['GenerateCombOutputFile'])
    controller.wf_GenerateStatFiles.set(data['ModelParams']['GenerateStatFiles'])

    wf_strata_tab.deconvert_tab_vars()




