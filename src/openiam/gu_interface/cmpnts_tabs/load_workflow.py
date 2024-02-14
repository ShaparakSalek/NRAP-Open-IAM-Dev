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
    aalf_tab, daa_tab, daaml_tab, fgaq_tab, fgaz_tab,
    ga_tab, atm_tab, psa_tab, strata_tab, cws_tab, AOR_wf_tab, locations)
from openiam.gu_interface.cmpnts_tabs.parameter_entry import ParameterEntry

from openiam.components.iam_base_classes import IAM_DIR

# Save location of source folder in the top level folder
USER_DIR = os.sep.join([IAM_DIR, 'examples', 'user'])

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
        elif key.find('Aquifer') != -1:
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

    if data['Workflow']['type'] == 'AoR':

        for key in ckeys:
            # Set up Workflow tab

            # Load figure resolution
            if key == 'Workflow':
                # Set figure DPI
                componentVars[key]['Params']['FigureDPI'].set(data[key]['Options']['FigureDPI'])

                # Set critical pressure strategy
                if data[key]['Options']['CriticalPressureMPa'] == 'Calculated':
                    componentVars[key]['Params']['CriticalPressureMPa'].set(True)
                    if 'OpenWellbore' in ckeys:
                        componentVars['OpenWellbore1']['Controls']['critPressureApproach'].set(True)
                else:
                    componentVars[key]['Params']['CriticalPressureMPa'].set(False)
                    crit_press = float(data[key]['Options']['CriticalPressureMPa'])
                    componentVars[key]['Params']['CriticalPressureMPa_no_calc'].set(crit_press)
                    if 'OpenWellbore' in ckeys:
                        componentVars['OpenWellbore1']['Controls']['critPressureApproach'].set(True)
                        componentVars['OpenWellbore1']['Controls']['enforceCritPressure'].set(True)
                        componentVars['OpenWellbore1']['Parameters']['critPressure']['value'] = crit_press * (10 ** 6)

                if 'AnalyticalReservoir' in ckeys:
                    if 'OpenWellbore' in ckeys or 'MultisegmentedWellbore' in ckeys:
                        componentVars[key]['Params']['BrineDensity'].set(data['AnalyticalReservoir']
                                                                             ['Parameters']
                                                                             ['brineDensity']
                                                                             ['value'])

            # Change brine density of components if necessary
            if controller.component_list['reservoir'][1].get() == 'AnalyticalReservoir' and \
                controller.component_list['wellbore'][1].get().replace(" ", "") in ['OpenWellbore',
                                                                                    'MultisegmentedWellbore']:
                brine_density = data['AnalyticalReservoir1']['Parameters']['brineDensity']['value']
                componentVars['Workflow']['Params']['BrineDensity'].set(brine_density)
                componentVars[resvar]['Params']['brineDensity']['value'].set(brine_density)
                componentVars[wellvar]['Params']['brineDensity']['value'].set(brine_density)

            # Load injection location information
            grid_names = ['xmin', 'xmax', 'xsize', 'ymin', 'ymax', 'ysize']
            for n in grid_names:
                componentVars['Workflow']['Locations']['grid'][n].set(data[wellvar]['Locations']['grid'][n])

            componentVars[wellvar]['Params']['number'] = int(data[wellvar]['Locations']['grid']['xsize'] *
                                                             data[wellvar]['Locations']['grid']['ysize'])

            # Set value for plotting injection sites
            plotInjSites = data['Workflow']['Options']['PlotInjectionSites']
            componentVars['Workflow']['Params']['PlotInjectionSites'].set(plotInjSites)

            if resvar == 'LookupTableReservoir1':
                inj_names = ['InjectionCoordx', 'InjectionCoordy']
                coord_names = ['xCoordinates', 'yCoordinates']
                for n in range(len(coord_names)):
                    componentVars['Workflow'][coord_names[n]].set(data['Workflow']['Options'][inj_names[n]])
