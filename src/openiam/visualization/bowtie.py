# -*- coding: utf-8 -*-
"""

Last Modified: January, 2024

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
LRST (Battelle|Leidos) supporting NETL
"""
import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from openiam.matk.sampleset import mean, std, percentile

from openiam.components.iam_base_classes import IAM_DIR

from .bowtie_label_setup import (UNIT_DICT, METRIC_LABEL_DICT, 
                                 METRIC_LABEL_CONCISE_DICT)

RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}

INPUT_ERR_MSG = ''.join([
    'The {} input for the Bowtie figure {} was given as {}. This ', 
    'input {}{}.'])

INPUT_COLS = ['ComponentNameList', 'NumberOfLocations', 'Metrics', 'ThresholdValues', 
              'Label', 'ReadResults', 'FileDirectory', 'FileName', 'Analysis', 
              'Realizations', 'Units']

BOOL_INPUT_COLS = ['ReadResults']

REQUIRED_INPUT_COLS = ['ComponentNameList', 'Metrics', 'ThresholdValues', 'Label']

REQUIRED_INPUT_WARNING_MSG = ''.join([
    'The input {} for the Bowtie figure {} was not found in the {} input file ', 
    'specified ({}). This input is required for the Bowtie figure.'])

STR_COLS = ['Label', 'Labels']

NUMERIC_COLS = ['ThresholdValue', 'ThresholdValues']

UNIT_ERR_MSG = ''.join([
    'For the Bowtie figure {}, the metric {} was given was for the {} index of ', 
    'the {} input. This metric was not recognized as an output type with units ', 
    'of kg/s. Check your input.'])

NEGATIVE_MASS_WARNING = ''.join([
    'While making the Bowtie figure {}, the {}{}{} included negative value(s) ', 
    '({}, in kg). Check your input.'])

EXPERT_PANEL = ['Expert Panel', 'ExpertPanel']

EXPERT_PANEL_VALS = ['Extremely Unlikely', 'Unlikely', 'Inconclusive', 
                     'Likely', 'Extremely Likely']

EXPERT_PANEL_WARNING = ''.join([
    'There was potentially an issue while assessing the {} input for the Bowtie plot {}. ', 
    'A value given for the {} index of the Metrics field was {}, which was ', 
    'not recognized as an expected value for Expert Panel feedback. The ', 
    'expected values are {}.'])

ANALYSIS_TYPE_LABEL = {'forward': 'Deterministic Simulation', 
                       'lhs': 'Stochastic Simulation\nwith {} Realizations', 
                       'parstudy': 'Stochastic Simulation\nwith {} Realizations'}

REALIZATION_ERROR_MSG = ''.join([
    'While making the Bowtie plot {}, there was an error. The analysis type ', 
    'was {}, but the ModelParams section of the .yaml file did not correctly ', 
    'specify the number of realizations used in the stochastic simulation. The ', 
    'number of realizations should be specified under "ModelParams: Analysis: siz".'])

BOWTIE_LEN_WARNING_MSG = ''.join([
    'The collection of results for the Bowtie figure {} may have encountered an ', 
    'error. During the collection of results for the {} of the figure, ', 
    'the number of collected result types was {} while it should have been {} ', 
    '(where the "correct" number is based on the {} input for the Bowtie figure). ', 
    'This situation can occur, for example, if a metric name or component name ', 
    'is misspelled. If a wellbore component is set to produce the metric ',
    'brine_aquifer1 but the Bowtie figure is set up to use brine_aquifer ', 
    '(lacking the number 1), then this error can occur. Otherwise, the component ', 
    'name could be entered as MultisegmentedWellbore1, when it was actually ', 
    'named msw1. Check your input.'])

INPUT_LEN_WARNING_MSG = ''.join([
    'There was an issue while handling the {} input for the Bowtie figure {}. ', 
    'The {} input was given{}, but this input did not have the same number of ', 
    'entries as the {} {}. This setup will likely cause an error. Check your ', 
    'input.'])

DATA_LOADING_ERR_MSG = ''.join([
    'While making the Bowtie plot {}, there was an error while attempting to ', 
    'read results from the file {}. The code was searching for data in a column ', 
    'labeled {}. The column was expected to have {} rows{}. Check your input.'])

TIME_LIST_ERR_MSG = ''.join([
    'The TimeList entry provided ({}, type: {}) for the Bowtie figure {} was not ', 
    'one of the expected input types. TimeList is expected to be "All" (for all ', 
    'times) or a list of times in years (e.g., [5, 10, 15, 20]). Check your ', 
    'input. The final model time will be used.'])

def bowtie_plot(yaml_data, model_data, sm, s, output_list, locations, output_dir, 
                name='Bowtie_Figure1', analysis='lhs', figsize=(12, 8), figdpi=100, 
                genfontsize=12, axislabelfontsize=14, titlefontsize=14, 
                boldlabels=True):
    """
    Makes a bowtie figure summarizing NRAP-Open-IAM results.

    :param yaml_data: Dictionary of input values
    :type yaml_data: dict

    :param model_data: Input from the 'ModelParams' section of the .yaml file
    :type model_data: dict

    :param locations: dictionary of locations assigned to each applicable component
    :type locations: dict

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param s: SampleSet with results of analysis performed on the system model.
        None for forward analysis.
    :type s: matk.SampleSet object

    :param name: Figure Name to be used/created.
    :type name: str

    :param analysis: Type of OpenIAM system analysis performed ('forward',
        'lhs', or 'parstudy')
    :type analysis: str

    :param savefig: Filename to save resulting figure to. No file saved if None.
    :type savefig: str

    :param figsize: width and height of the figure (width, height), in inches.
        Default value is (10, 8).
    :type figsize: tuple

    :param genfontsize: fontsize for tick labels, etc.
    :type genfontsize: float or int

    :param axislabelfontsize: fontsize for x and y axis labels
    :type axislabelfontsize: float or int

    :param titlefontsize: fontsize for the title
    :type titlefontsize: float or int

    :param boldlabels: option to use bold x and y labels and bold titles. Set to
        True for bold labels, False for normal labels.
    :type boldlabels: bool

    :return: None
    """
    yaml_input = yaml_data['Plots'][name]['Bowtie']
    
    contribs = yaml_input['Contributors']
    conseqs = yaml_input['Consequences']
    
    top_event = yaml_input.get('TopEventLabel', None)
    
    selected_times = yaml_input.get('TimeList', None)
    
    if not selected_times is None and selected_times != 'All':
        selected_times = np.array(selected_times)
    
    figdpi = yaml_input.get('FigureDPI', figdpi)

    save_results = yaml_input.get('SaveCSVFiles', True)
    
    contribs = check_bowtie_input(name, contribs, 'Contributors')
    
    conseqs = check_bowtie_input(name, conseqs, 'Consequences')
    
    if analysis == 'forward':
        realizations = None
    elif analysis in ['lhs', 'parstudy']:
        realizations = model_data['Analysis'].get('siz', None)
        
        if realizations is None:
            err_msg = REALIZATION_ERROR_MSG.format(name, analysis)
            
            logging.error(err_msg)
            raise KeyError(err_msg)
        else:
            realizations = int(realizations)
    
    results_contribs = get_bowtie_results(
        sm, s, output_list, contribs, 'Contributors', name, analysis=analysis, 
        output_dir=output_dir, realizations=realizations)
    
    results_conseqs = get_bowtie_results(
        sm, s, output_list, conseqs, 'Consequences', name, analysis=analysis, 
        output_dir=output_dir, realizations=realizations)
    
    time_array = sm.time_array / 365.25
    
    make_bowtie_figure(name, output_dir, contribs, conseqs, results_contribs, 
                       results_conseqs, top_event,  time_array, 
                       selected_times=selected_times, figsize=figsize, 
                       figdpi=figdpi, boldlabels=boldlabels, 
                       text_font_size=genfontsize, 
                       text_font_size_2=(genfontsize * 0.75), 
                       title_font_size=titlefontsize, save_results=save_results, 
                       analysis=analysis, realizations=realizations)


def check_bowtie_input(name, bowtie_input, input_type):
    """
    Checks the input provided for 'Contributors' or 'Consequences'. Returns a 
    dictionary containing the input.
    """
    if isinstance(bowtie_input, (str, dict)):
        if isinstance(bowtie_input, str):
            file_path = os.path.join(IAM_DIR, bowtie_input)
            if os.path.isfile(file_path):
                data = pd.read_csv(
                    file_path, delimiter=",", dtype='str', comment='#')
            else:
                logging.error(INPUT_ERR_MSG.format(
                    input_type, name, file_path, 'file ', 'was not found'))
        
        input_dict = dict()
        for key in INPUT_COLS:
            if key in data:
                if not key in input_dict:
                    input_dict[key] = []
                
                for ind, val in enumerate(data[key]):
                    input_dict[key].append([])
                    
                    if isinstance(val, str):
                        # The input for the keys in STR_COLS will always be string 
                        # input that is not intended to be converted into a list.
                        # For example, the input for 'Label' is always a string. 
                        if not key in STR_COLS:
                            if ',' in val:
                                edited_vals = val.split(',')
                                
                                for ind, edited_val in enumerate(edited_vals):
                                    # Get rid of blank spaces
                                    edited_vals[ind] = evaluate_input_val(
                                        key, edited_val.replace(' ', '')) 
                                
                                input_dict[key][-1] = edited_vals
                                
                            else:
                                input_dict[key][-1].append(evaluate_input_val(key, val))
                        else:
                            input_dict[key][-1].append(evaluate_input_val(key, val))
                    else:
                        input_dict[key][-1].append(evaluate_input_val(key, val))
                    
            elif (not key in data) and (key in REQUIRED_INPUT_COLS):
                err_msg = REQUIRED_INPUT_WARNING_MSG.format(
                    key, name, input_type, file_path)
                logging.error(err_msg)
                
    else:
        err_msg = INPUT_ERR_MSG.format(
            input_type, name, bowtie_input, '', 
            'was not recognized as a suitable input format (.csv file or dictionary)')
        
        logging.error(err_msg)
        
        raise KeyError(err_msg)
    
    return input_dict


def evaluate_input_val(key, value):
    """
    Checks if the key is in the list NUMERIC_COLS. The keys in that list are 
    generally meant to be taken as numeric values rather than strings.
    """
    if key in NUMERIC_COLS:
        try:
            return eval(value)
        except:
            return None
    elif value in ['None', 'none']:
        return None
    elif value in ['False', 'FALSE', 'false']:
        return False
    elif value in ['True', 'TRUE', 'true']:
        return True
    else:
        if isinstance(value, str):
            return value
        elif np.isnan(value):
            return None
        else:
            return value


def get_optional_bowtie_input(name, bowtie_input, metric_list, bowtie_input_type, 
                              optional_input_type):
    """
    Checks the input file for a given entry type. If present, those inputs are 
    returned. Otherwise, the funciton returns a list of None values.
    """
    input_list = [None] * len(metric_list)
    
    if optional_input_type in bowtie_input:
        input_list = bowtie_input[optional_input_type]
        
        if len(input_list) != len(metric_list):
            err_msg = INPUT_LEN_WARNING_MSG.format(
                bowtie_input_type, name, optional_input_type, '', 
                'Metrics', 'input')
            
            logging.error(err_msg)
            
            input_list = [False] * len(metric_list)
        else:
            for row_ind in range(len(input_list)):
                if isinstance(input_list[row_ind], list):
                    for metric_ind, input_val in enumerate(input_list[row_ind]):
                        if optional_input_type in BOOL_INPUT_COLS:
                            input_list[row_ind][metric_ind] = eval(
                                input_list[row_ind][metric_ind])
                        else:
                            input_list[row_ind][metric_ind] = input_list[
                                row_ind][metric_ind]
                else:
                    if optional_input_type in BOOL_INPUT_COLS:
                        input_list[row_ind] = eval(input_list[row_ind])
                    else:
                        input_list[row_ind] = input_list[row_ind]
    
    return input_list


def read_results_file(time_array, full_obs_nm, file_dir, file_name, bowtie_input, 
                      input_type, name, realizations=None, analysis='forward', 
                      reals=None, analysis_type=None):
    """
    Function that reads a .csv file and returns the results for a particular 
    metric. If the simulation is deterministic (analysis = 'forward'), then 
    the results should have the same length as the time_array. When using a 
    deterministic model, all returned variables will be None except for values.
    Note that the 'reals' and 'analysis_type' variables are taken from input 
    provided through the 'Contributors' or 'Consequences' input files. These 
    variables are meant to be used for output from a different program. For 
    example, the NRAP-Open-IAM simulation is stochastic (analysis = 'lhs' or 
    analysis = 'parstudy'), but the results taken from another program may not 
    be able to be produced with the same approach (e.g., only deterministic). 
    These variables are provided to allow for greater flexibility when reading 
    values into the Bowtie plot.
    """
    file_path = os.path.join(file_dir, file_name)
    
    if os.path.isfile(file_path):
        data = pd.read_csv(
            file_path, delimiter=",", dtype='str', comment='#')
    else:
        logging.error(INPUT_ERR_MSG.format(
            input_type, name, file_path, 'file ', 'was not found'))
    
    # The variable 'realizations' refers to the realizations used in NRAP-Open-IAM. 
    # The user might want to use output from a different model (e.g., a cost model), 
    # and that model might have a different setup. The variables 'analysis_type' 
    # and 'reals' refer to the approach used by a separate model. For example, 
    # the NRAP-Open-IAM results might use LHS while the separate model can only 
    # use a deterministic analysis type (forward). This plot type is meant to 
    # be flexibile for the user.
    if not analysis_type is None:
        analysis = analysis_type
    
    if not reals is None:
        realizations = int(reals)
    
    if analysis == 'forward':
        try:
            values = data[full_obs_nm].values
        except:
            err_msg = DATA_LOADING_ERR_MSG.format(
                name, file_path, full_obs_nm, len(time_array), 
                ' to align with the time range used')
            
            logging.error(err_msg)
            raise KeyError(err_msg)
        
        values = np.array([float(val) for val in values])
        values = values[~np.isnan(values)]
        
        mean_vals = None
        std_vals = None
        min_vals = None
        max_vals = None
        
        if len(values) != len(time_array):
            err_msg = INPUT_LEN_WARNING_MSG.format(
                input_type, name, full_obs_nm, ' as being stored in the file {}'.format(
                    file_path), 'time range', 
                'used (values were: {}; had a length of {},'.format(
                    values, len(values))
                + ' should have a length of {} to match the time range: {})'.format(
                    len(time_array), time_array))
            
            logging.error(err_msg)
            raise KeyError(err_msg)
        
    elif analysis in ['lhs', 'parstudy']:
        values = np.zeros((len(time_array), realizations))
        
        for time_index in range(1, values.shape[0]):
            for real in range(realizations):
                col_name = 'Realization {}'.format(real)
                
                try:
                    vals_temp = data[col_name].values
                except:
                    err_msg = DATA_LOADING_ERR_MSG.format(
                        name, file_path, col_name, len(time_array), 
                        ' to align with the time range used')
                    
                    logging.error(err_msg)
                    raise KeyError(err_msg)
                
                vals_temp = np.array([float(val) for val in vals_temp])
                vals_temp = vals_temp[~np.isnan(vals_temp)]
                
                if len(vals_temp) != len(time_array):
                    err_msg = INPUT_LEN_WARNING_MSG.format(
                        input_type, name, col_name + ' for the ComponentName.Metric ' 
                        + full_obs_nm + ' values', ' as being stored in the file {}'.format(
                            file_path), 'time range', 
                        'used (values were: {}; had a length of {},'.format(
                            vals_temp, len(vals_temp)) 
                        + ' should have a length of {} to match the time range: {})'.format(
                            len(time_array), time_array))
                    
                    logging.error(err_msg)
                    raise KeyError(err_msg)
                
                values[:, real] = vals_temp
                
                del vals_temp
        
        mean_vals = np.zeros(len(time_array))
        min_vals = np.zeros(len(time_array))
        max_vals = np.zeros(len(time_array))
        std_vals = np.zeros(len(time_array))
        
        for time_index in range(1, len(values)):
            mean_vals[time_index] = mean(values[time_index,:])
            min_vals[time_index] = np.min(values[time_index,:])
            max_vals[time_index] = np.max(values[time_index,:])
            std_vals[time_index] = std(values[time_index,:])
    
    return values, mean_vals, min_vals, max_vals, std_vals


def get_current_entry(input_from_file, row_ind, metric_ind):
    """
    Checks the lists produced by get_optional_bowtie_input() and returns the 
    appropriate values, depending on the indices row_ind and metric_ind.
    """
    if isinstance(input_from_file[row_ind], list):
        if len(input_from_file[row_ind]) == 1:
            current_input = input_from_file[row_ind][0]
        else:
            current_input = input_from_file[row_ind][metric_ind]
    else:
        current_input = input_from_file[row_ind]
    
    return current_input


def get_bowtie_results(sm, s, output_list, bowtie_input, input_type, name, 
                       analysis='forward', output_dir=None, realizations=None):
    """
    Compiles the results required for the Bowtie figure. Results can be extracted 
    from the NRAP-Open-IAM system model (using sm, s, and output_list) or read 
    from pre-made files.
    """
    time_array = sm.time_array / 365.25
    
    component_list = bowtie_input['ComponentNameList']
    metric_list = bowtie_input['Metrics']
    
    # The following lists are only used for results read from pre-made files
    read_results_list = get_optional_bowtie_input(
        name, bowtie_input, metric_list, input_type, 'ReadResults')
    
    file_dir_list = get_optional_bowtie_input(
        name, bowtie_input, metric_list, input_type, 'FileDirectory')
    
    file_name_list = get_optional_bowtie_input(
        name, bowtie_input, metric_list, input_type, 'FileName')
    
    analysis_type_list = get_optional_bowtie_input(
        name, bowtie_input, metric_list, input_type, 'Analysis')
    
    realizations_list = get_optional_bowtie_input(
        name, bowtie_input, metric_list, input_type, 'Realizations')
    
    units_list = get_optional_bowtie_input(
        name, bowtie_input, metric_list, input_type, 'Units')
    
    results = []
    
    for row_ind, (comps, metrics) in enumerate(zip(component_list, metric_list)):
        results.append({})
        
        for metric_ind, (comp, output_nm) in enumerate(zip(comps, metrics)):
            
            # Check whether the results are read from a file
            read_check = get_current_entry(read_results_list, row_ind, metric_ind)
            
            if read_check:
                # Get the file directory, if one is used for reading files
                file_dir_temp = get_current_entry(file_dir_list, row_ind, metric_ind)
                
                # Get the file directory, if one is used for reading files
                if file_dir_temp:
                    file_dir = os.path.join(IAM_DIR, file_dir_temp)
                else:
                    file_dir = os.path.join(IAM_DIR, output_dir, 'csv_files')
                
                # Get the file name, if one is used for reading files
                file_name = get_current_entry(file_name_list, row_ind, metric_ind)
                
                # Get the inputs 'Analysis', 'Realizations', and 'Units'. These 
                # are meant for metrics that are not produced by NRAP-Open-IAM.
                analysis_type = get_current_entry(analysis_type_list, 
                                                  row_ind, metric_ind)
                reals = get_current_entry(realizations_list, 
                                          row_ind, metric_ind)
                units = get_current_entry(units_list, row_ind, metric_ind)
                
                # If the metric type has not been added as a key, add it
                if not output_nm in results[row_ind]:
                    results[row_ind][output_nm] = {}
                
                # Dictionary for the metric type, for the specific component 
                # being handled
                results[row_ind][output_nm][comp] = {}
                
                results = check_results(
                    sm, s, results, row_ind, metric_ind, comp, 
                    output_nm, time_array, bowtie_input, input_type, name, 
                    analysis=analysis, read_results=read_check, 
                    file_dir=file_dir, file_name=file_name, 
                    realizations=realizations, analysis_type=analysis_type, 
                    reals=reals, units=units)
                
            else:
                if comp in EXPERT_PANEL:
                    if not output_nm in EXPERT_PANEL_VALS:
                        warning_msg = EXPERT_PANEL_WARNING.format(
                            input_type, name, row_ind, output_nm, EXPERT_PANEL_VALS)
                        
                        logging.warning(warning_msg)

                    results[row_ind][EXPERT_PANEL[0]] = output_nm
                    
                else:
                    for output_component in list(output_list.keys()):
                        # Checks if the output is produced by this component and if the 
                        # component name is in the list given for Contributor / Consequence.
                        applies = (output_nm in output_list[output_component] 
                                   and comp in output_component.name)
                        
                        if applies:
                            # If the metric type has not been added as a key, add it
                            if not output_nm in results[row_ind]:
                                results[row_ind][output_nm] = {}
                            
                            # Dictionary for the metric type, for the specific component 
                            # being handled
                            results[row_ind][output_nm][output_component.name] = {}
                            
                            results = check_results(
                                sm, s, results, row_ind, metric_ind, output_component, 
                                output_nm, time_array, bowtie_input, input_type, name, 
                                analysis=analysis, realizations=realizations)
    
    return results


def check_results(sm, s, results, row_ind, metric_ind, output_component, 
                  output_nm, time_array, bowtie_input, input_type, name, 
                  analysis='forward', read_results=False, file_dir=None, 
                  file_name=None, realizations=None, analysis_type=None, 
                  reals=None, units=None):
    """
    """
    thresh_vals = bowtie_input['ThresholdValues'][row_ind]
    
    if len(thresh_vals) > 1:
        thresh = thresh_vals[metric_ind]
    else:
        thresh = thresh_vals[0]
    
    if read_results:
        full_obs_nm = output_component + '.' + output_nm
        
        values, mean_vals, min_vals, max_vals, std_vals = read_results_file(
            time_array, full_obs_nm, file_dir, file_name, bowtie_input, 
            input_type, name, realizations=realizations, analysis=analysis, 
            reals=reals, analysis_type=analysis_type)
        
        prob_vals = get_prob_vals(values, thresh, analysis=analysis, 
                                  analysis_type=analysis_type)
    else:
        values, mean_vals, std_vals, min_vals, max_vals, prob_vals = get_results(
            sm, s, time_array, output_component, output_nm, thresh, 
            row_ind, name, input_type, analysis=analysis)
        
        # These are only needed for results read from pre-made files
        analysis_type = None
        reals = None
        units = None
    
    if isinstance(output_component, str):
        comp_name = output_component
    else:
        comp_name = output_component.name
    
    results[row_ind][output_nm][comp_name] = {
        'values': values, 'mean_vals': mean_vals, 'std_vals': std_vals, 
        'min_vals': min_vals, 'max_vals': max_vals, 'prob_vals': prob_vals, 
        'thresh': thresh, 'analysis': analysis_type, 'reals': reals, 
        'units': units}
    
    return results


def get_results(sm, s, time_array, output_component, output_nm, thresh, row_ind, 
                name, input_type, analysis='forward'):
    """
    Extracts results for an NRAP-Open-IAM simulation using sm and s.
    """
    if analysis == 'forward':
        values = sm.collect_observations_as_time_series(
            output_component, output_nm)
        
        mean_vals = None
        std_vals = None
        min_vals = None
        max_vals = None
        
    elif analysis in ['lhs', 'parstudy']:
        ind_list = list(range(len(time_array)))
        
        full_obs_nm = '.'.join([output_component.name, output_nm])
        obs_names = [full_obs_nm + '_{0}'.format(indd)
                     for indd in ind_list]
        
        values = np.array(
            [s.recarray[full_obs_nm + '_' + str(indd)] for indd in ind_list])
        
        mean_vals = mean(s.recarray[obs_names])
        std_vals = std(s.recarray[obs_names])
        
        obs_percentiles = percentile(s.recarray[obs_names], [0, 100])
        
        min_vals = obs_percentiles[0, :]
        max_vals = obs_percentiles[1, :]
    
    prob_vals = get_prob_vals(values, thresh, analysis=analysis)
    
    return values, mean_vals, std_vals, min_vals, max_vals, prob_vals


def get_prob_vals(values, thresh, analysis='forward', analysis_type=None):
    """
    Calculates the probability of whether a given threshold is exceeded in a 
    value array. The probability is taken over time, with later times reflecting 
    if the threshold was exceeded at an earlier time. If the analysis type is 
    forward (deterministic), the values array has a length equal to the number 
    of time steps. If the analysis type is lhs or parstudy (stochastic), then 
    the values.shape is (number of time steps, number of realizations).
    """
    # The variable 'analysis_type' refers to the approach used by a separate
    # model. A separate model can use a different approach than the one used 
    # by NRAP-Open-IAM (which is described by the analysis variable).
    if not analysis_type is None:
        analysis = analysis_type
    
    if analysis == 'forward':
        prob_vals = values.copy() * 0
        
        for t_index in range(1, len(values)):
            if np.max(values[:t_index]) >= thresh:
                prob_vals[t_index] = 1
        
    elif analysis in ['lhs', 'parstudy']:
        prob_vals = np.zeros(values.shape[0])
        
        # Get the probability, taken as the number of realizations in 
        # which the threshold is exceeded divided by the total number 
        # of realizations.
        for t_index in range(1, values.shape[0]):
            for real in range(values.shape[1]):
                if np.max(values[:t_index, real]) >= thresh:
                    prob_vals[t_index] += 1
        
        prob_vals /= values.shape[1]
    
    prob_vals_edit = np.zeros(values.shape[0])
    
    # The probability at each time should reflect the preceding times
    for t_index in range(1, len(prob_vals)):
        prob_vals_edit[t_index] = np.max(prob_vals[:t_index])
    
    prob_vals = prob_vals_edit.copy()
    
    return prob_vals


def make_bowtie_figure(name, output_dir, contribs, conseqs, results_contribs, 
                       results_conseqs, top_event, time_array, selected_times=None, 
                       plot_metric_labels=False, show_std_option=True, 
                       plot_legend=True, plot_thresholds=True, figsize=(15, 12), 
                       figdpi=100, fontname='Arial', gen_font_size=10, 
                       text_font_size=12, text_font_size_2=9, title_font_size=12, 
                       colormap=None, boldlabels=True, save_results=True, 
                       analysis='forward', realizations=None):
    """
    Produces a bowtie figure.
    """
    if not selected_times is None:
        if isinstance(selected_times, (list, np.ndarray)):
            time_indices = []
            
            for time in selected_times:
                time_abs_diff = np.abs(time_array - time)
                closest_time_index = np.argmin(time_abs_diff)
                time_indices.append(closest_time_index)
            
        elif isinstance(selected_times, str):
            if selected_times == 'All':
                selected_times = np.arange(0, len(time_array)).tolist()
            else:
                err_msg = TIME_LIST_ERR_MSG.format(
                    selected_times, type(selected_times), name)
                
                logging.error(err_msg)
                
                time_indices = [len(time_array) - 1]
        else:
            err_msg = TIME_LIST_ERR_MSG.format(
                selected_times, type(selected_times), name)
            
            logging.error(err_msg)
            
            time_indices = [len(time_array) - 1]
    else:
        time_indices = [len(time_array) - 1]
    
    if analysis in ['lhs', 'parstudy'] and realizations is None:
        err_msg = REALIZATION_ERROR_MSG.format(name, analysis)
        
        logging.error(err_msg)
        raise KeyError(err_msg)
    
    # Figures
    font = RC_FONT
    font['size'] = gen_font_size
    plt.rc('font', **font)
    
    if boldlabels:
        fweight = 'bold'
    else:
        fweight = 'normal'
    
    text_zorder = 2.0e+6
    
    label_fcolor = 'k'
    
    whisker_color = 'k'
    lwidth = 2
    lstyle = '-'
    whisker_zorder = 2
    
    # IAM_bbox = {}
    # Expert_bbox = {}
    
    # center position of the bowtie figure
    center_x = 0
    center_y = 0
    
    circle_radius = 2
    circ_edge_color = 'k'
    circ_face_color = 'firebrick'
    circ_lwidth = 2
    circle_zorder = 3
    
    top_event_line_len = 12.5
    top_event_text_ydiff = top_event_line_len * 0.05
    top_event_whisker_color = 'firebrick'
    
    # y limits for the whiskers
    whisker_min_y = -10
    whisker_max_y = 10
    
    # x length of the whisker, where it remains flat (not angled)
    whisker_len_x_straight = 15
    
    # x length of the region where whiskers angle towards the center position
    whisker_len_x_angled = 5
    
    contribs_x = np.array([
        center_x - (whisker_len_x_straight + whisker_len_x_angled), 
        center_x - ((whisker_len_x_straight * 0.25) + whisker_len_x_angled), 
        center_x - whisker_len_x_angled])
    
    conseqs_x = np.array([
        center_x + whisker_len_x_angled, 
        center_x + ((whisker_len_x_straight *  0.25) + whisker_len_x_angled), 
        center_x + (whisker_len_x_straight + whisker_len_x_angled)])
    
    contribs_rows = len(results_contribs)
    conseqs_rows = len(results_conseqs)
    
    contribs_y = np.linspace(whisker_min_y, whisker_max_y, num=contribs_rows)
    conseqs_y = np.linspace(whisker_min_y, whisker_max_y, num=conseqs_rows)
    
    # y space used to offset the text from the whisker line
    contribs_text_ydiff = np.abs(contribs_y[1] - contribs_y[0]) * 0.1
    conseqs_text_ydiff = np.abs(conseqs_y[1] - conseqs_y[0]) * 0.1
    
    # Flip these so that the order you enter risks into the file is the order
    # order they are displayed on the figure.
    contribs_y = np.flip(contribs_y)
    conseqs_y = np.flip(conseqs_y)
    
    time_x_offset = whisker_len_x_angled + whisker_len_x_straight
    time_y_offset = 15
    
    analysis_x_offset = (whisker_len_x_angled + (whisker_len_x_straight * 0.25))
    analysis_y_offset = 15
    
    # Position of the legend
    legend_x_offset = (whisker_len_x_straight * 0.5)
    legend_line_x = [center_x + whisker_len_x_angled, 
                     center_x + whisker_len_x_angled + (legend_x_offset * 0.5), 
                     center_x + whisker_len_x_angled + legend_x_offset]
    
    legend_y_offset = 15
    legend_line_y = [center_y - legend_y_offset, 
                     center_y - legend_y_offset, center_y - legend_y_offset]
    
    prob_legend_y_diff = 0.67
    val_legend_y_diff = 0.33
    
    thresh_legend_x = center_x + (whisker_len_x_straight + whisker_len_x_angled)
    thresh_legend_y = center_y - legend_y_offset
    
    thresh_text_color = 'firebrick'
    
    expert_panel_color = 'goldenrod'
    
    expert_legend_x = thresh_legend_x
    expert_legend_y = thresh_legend_y - 0.75
    
    for t_index in time_indices:
        # This indicates if any of the metrics used had ranges from multiple 
        # locations. If multiple locations are used but they all had the same 
        # value for the metric, this will remain as False.
        mult_loc_check_overall = False
        
        # Get the values stored in the results dictionaries
        contribs_metrics = []
        contribs_vals = []
        contribs_min_vals = []
        contribs_max_vals = []
        contribs_std_vals = []
        contribs_probs = []
        contribs_threshs = []
        contribs_analysis = []
        contribs_reals = []
        contribs_units = []
        
        for row in range(len(results_contribs)):
            max_prob = 0
            
            for key1 in results_contribs[row].keys():
                contribs_metrics.append(key1)
                
                if key1 is EXPERT_PANEL[0]:
                    max_prob = results_contribs[row][key1]
                    vals = None
                    min_vals = None
                    max_vals = None
                    std_vals = None
                    thresh = None
                else:
                    vals = []
                    min_vals = []
                    max_vals = []
                    std_vals = []
                    for key2 in results_contribs[row][key1].keys():
                        max_prob = np.max([max_prob, 
                                           results_contribs[
                                               row][key1][key2]['prob_vals'][t_index]])
                        
                        thresh = results_contribs[row][key1][key2]['thresh']
                        analysis_type = results_contribs[row][key1][key2]['analysis']
                        reals = results_contribs[row][key1][key2]['reals']
                        unit = results_contribs[row][key1][key2]['units']
                        
                        if not analysis_type is None:
                            analysis_temp = analysis_type
                        else:
                            analysis_temp = analysis
                        
                        if analysis_temp == 'forward':
                            vals.append(results_contribs[
                                row][key1][key2]['values'][t_index])
                            min_vals.append(None)
                            max_vals.append(None)
                            std_vals.append(None)
                        if analysis_temp in ['lhs', 'parstudy']:
                            vals.append(results_contribs[
                                row][key1][key2]['mean_vals'][t_index])
                            min_vals.append(results_contribs[
                                row][key1][key2]['min_vals'][t_index])
                            max_vals.append(results_contribs[
                                row][key1][key2]['max_vals'][t_index])
                            std_vals.append(results_contribs[
                                row][key1][key2]['std_vals'][t_index])
                
                contribs_probs.append(max_prob)
                contribs_vals.append(vals)
                contribs_min_vals.append(min_vals)
                contribs_max_vals.append(max_vals)
                contribs_std_vals.append(std_vals)
                contribs_threshs.append(thresh)
                contribs_analysis.append(analysis_type)
                contribs_reals.append(reals)
                contribs_units.append(unit)
        
        conseqs_metrics = []
        conseqs_vals = []
        conseqs_min_vals = []
        conseqs_max_vals = []
        conseqs_std_vals = []
        conseqs_probs = []
        conseqs_threshs = []
        conseqs_analysis = []
        conseqs_reals = []
        conseqs_units = []
        
        for row in range(len(results_conseqs)):
            max_prob = 0
            
            for key1 in results_conseqs[row].keys():
                conseqs_metrics.append(key1)
                
                if key1 is EXPERT_PANEL[0]:
                    max_prob = results_conseqs[row][key1]
                    vals = None
                    min_vals = None
                    max_vals = None
                    std_vals = None
                    thresh = None
                else:
                    vals = []
                    min_vals = []
                    max_vals = []
                    std_vals = []
                    for key2 in results_conseqs[row][key1].keys():
                        max_prob = np.max([max_prob, 
                                           results_conseqs[
                                               row][key1][key2]['prob_vals'][t_index]])
                        
                        thresh = results_conseqs[row][key1][key2]['thresh']
                        analysis_type = results_conseqs[row][key1][key2]['analysis']
                        reals = results_conseqs[row][key1][key2]['reals']
                        unit = results_conseqs[row][key1][key2]['units']
                        
                        if not analysis_type is None:
                            analysis_temp = analysis_type
                        else:
                            analysis_temp = analysis
                        
                        if analysis_temp == 'forward':
                            vals.append(results_conseqs[
                                row][key1][key2]['values'][t_index])
                            min_vals.append(None)
                            max_vals.append(None)
                            std_vals.append(None)
                        if analysis_temp in ['lhs', 'parstudy']:
                            vals.append(results_conseqs[
                                row][key1][key2]['mean_vals'][t_index])
                            min_vals.append(results_conseqs[
                                row][key1][key2]['min_vals'][t_index])
                            max_vals.append(results_conseqs[
                                row][key1][key2]['max_vals'][t_index])
                            std_vals.append(results_conseqs[
                                row][key1][key2]['std_vals'][t_index])
                
                conseqs_probs.append(max_prob)
                conseqs_vals.append(vals)
                conseqs_min_vals.append(min_vals)
                conseqs_max_vals.append(max_vals)
                conseqs_std_vals.append(std_vals)
                conseqs_threshs.append(thresh)
                conseqs_analysis.append(analysis_type)
                conseqs_reals.append(reals)
                conseqs_units.append(unit)
        
        if len(contribs_metrics) != len(results_contribs):
            err_msg = BOWTIE_LEN_WARNING_MSG.format(
                name, 'Contributors', len(contribs_metrics), len(results_contribs), 
                'Contributors')
            logging.error(err_msg)
        
        if len(conseqs_metrics) != len(results_conseqs):
            err_msg = BOWTIE_LEN_WARNING_MSG.format(
                name, 'Consequences', len(conseqs_metrics), len(results_conseqs), 
                'Consequences')
            logging.error(err_msg)
        
        fig = plt.figure(figsize=figsize, dpi=figdpi)

        ax = plt.gca()
        
        circle = plt.Circle((center_x, center_y), circle_radius, 
                            linewidth=circ_lwidth, edgecolor=circ_edge_color, 
                            color=circ_face_color, zorder=circle_zorder)
        
        ax.add_patch(circle)
        
        ax.text(center_x, center_y, 'Top Event', fontsize=text_font_size, 
                fontweight=fweight, color=label_fcolor, zorder=text_zorder, 
                horizontalalignment='center', verticalalignment='center')
        
        if top_event:
            ax.plot([center_x, center_x], [center_y, center_y + top_event_line_len], 
                    color=top_event_whisker_color, zorder=whisker_zorder)
            
            ax.text(
                center_x, center_y + top_event_line_len + top_event_text_ydiff, 
                top_event, fontsize=text_font_size, fontweight=fweight, 
                color=label_fcolor, zorder=text_zorder, 
                horizontalalignment='center', verticalalignment='center')
        
        # Plot Contributors
        for row, y in enumerate(contribs_y):
            ax.plot([contribs_x[0], contribs_x[2]], [y, y], color=whisker_color, 
                    linewidth=lwidth, linestyle=lstyle, zorder=whisker_zorder)
            
            ax.plot([contribs_x[2], center_x], [y, center_y], color=whisker_color, 
                    linewidth=lwidth, linestyle=lstyle, zorder=whisker_zorder)
            
            ax.text(contribs_x[0], y + contribs_text_ydiff, 
                    contribs['Label'][row][0], fontsize=text_font_size, 
                    fontweight=fweight, zorder=text_zorder, 
                    horizontalalignment='center', verticalalignment='center')
            
            if plot_metric_labels:
                if contribs_metrics[row] in METRIC_LABEL_DICT:
                    metric_label = METRIC_LABEL_DICT[contribs_metrics[row]]
                else:
                    metric_label = contribs_metrics[row]
                
                ax.text(contribs_x[0], y - (conseqs_text_ydiff * 1.25), metric_label, 
                        fontsize=text_font_size_2, fontweight=fweight, zorder=text_zorder, 
                        horizontalalignment='center', verticalalignment='center')
            
            if plot_thresholds:
                y_thresh_offset = 0
                if plot_metric_labels:
                    # This is menat to prevent the text from overlapping
                    y_thresh_offset += contribs_text_ydiff
                
                if not (contribs_metrics[row] in EXPERT_PANEL):
                    if contribs_metrics[row] in METRIC_LABEL_CONCISE_DICT:
                        metric_label = METRIC_LABEL_CONCISE_DICT[contribs_metrics[row]]
                    else:
                        metric_label = contribs_metrics[row]
                    
                    thresh_text = get_val_str(contribs_threshs[row], None, 
                                              contribs_metrics[row], 
                                              show_std_option=False, 
                                              analysis_type=contribs_analysis[row], 
                                              unit=contribs_units[row])
                    
                    thresh_label = metric_label + ' \u2265 ' + thresh_text
                    label_color = thresh_text_color
                else:
                    thresh_label = 'Expert Panel Review'
                    label_color = expert_panel_color
                
                ax.text(contribs_x[0], y - (contribs_text_ydiff * 1.25) - y_thresh_offset, 
                        thresh_label, color=label_color, 
                        fontsize=text_font_size_2, fontweight=fweight, 
                        zorder=text_zorder, horizontalalignment='center', 
                        verticalalignment='center')
            
            prob_text, val_text, mult_loc_check = get_result_text(
                contribs_probs[row], contribs_vals[row], contribs_std_vals[row], 
                contribs_metrics[row], show_std_option=show_std_option, 
                analysis_type=contribs_analysis[row], unit=contribs_units[row])
            
            if mult_loc_check:
                mult_loc_check_overall = mult_loc_check
                
            ax.text(contribs_x[1], y + contribs_text_ydiff, prob_text, 
                    fontsize=text_font_size, fontweight=fweight, zorder=text_zorder, 
                    horizontalalignment='center', verticalalignment='center')
            
            if val_text:
                ax.text(contribs_x[1], y - (contribs_text_ydiff * 0.9), val_text, 
                        fontsize=text_font_size_2, fontweight=fweight, zorder=text_zorder, 
                        horizontalalignment='center', verticalalignment='top')
        
        # Plot Consequences
        for row, y in enumerate(conseqs_y):
            ax.plot([conseqs_x[0], conseqs_x[2]], [y, y], color=whisker_color, 
                    linewidth=lwidth, linestyle=lstyle, zorder=whisker_zorder)
            
            ax.plot([conseqs_x[0], center_x], [y, center_y], color=whisker_color, 
                    linewidth=lwidth, linestyle=lstyle, zorder=whisker_zorder)
            
            ax.text(conseqs_x[2], y + conseqs_text_ydiff, conseqs['Label'][row][0], 
                    fontsize=text_font_size, fontweight=fweight, zorder=text_zorder, 
                    horizontalalignment='center', verticalalignment='center')
            
            if plot_metric_labels:
                if conseqs_metrics[row] in METRIC_LABEL_DICT:
                    metric_label = METRIC_LABEL_DICT[conseqs_metrics[row]]
                else:
                    metric_label = conseqs_metrics[row]
                
                ax.text(conseqs_x[2], y - (conseqs_text_ydiff * 1.25), metric_label, 
                        fontsize=text_font_size_2, fontweight=fweight, zorder=text_zorder, 
                        horizontalalignment='center', verticalalignment='center')
            
            if plot_thresholds:
                y_thresh_offset = 0
                if plot_metric_labels:
                    # This is menat to prevent the text from overlapping
                    y_thresh_offset += conseqs_text_ydiff
                
                if not (conseqs_metrics[row] in EXPERT_PANEL):
                    if conseqs_metrics[row] in METRIC_LABEL_CONCISE_DICT:
                        metric_label = METRIC_LABEL_CONCISE_DICT[conseqs_metrics[row]]
                    else:
                        metric_label = conseqs_metrics[row]
                    
                    thresh_text = get_val_str(conseqs_threshs[row], None, 
                                              conseqs_metrics[row], 
                                              show_std_option=False, 
                                              analysis_type=conseqs_analysis[row], 
                                              unit=conseqs_units[row])
                    
                    thresh_label = metric_label + ' \u2265 ' + thresh_text
                    label_color = thresh_text_color
                else:
                    thresh_label = 'Expert Panel Review'
                    label_color = expert_panel_color
                
                ax.text(conseqs_x[2], y - (conseqs_text_ydiff * 1.25) - y_thresh_offset, 
                        thresh_label, color=label_color, 
                        fontsize=text_font_size_2, fontweight=fweight, 
                        zorder=text_zorder, horizontalalignment='center', 
                        verticalalignment='center')
            
            prob_text, val_text, mult_loc_check = get_result_text(
                conseqs_probs[row], conseqs_vals[row], conseqs_std_vals[row], 
                conseqs_metrics[row], show_std_option=show_std_option, 
                analysis_type=conseqs_analysis[row], unit=conseqs_units[row])
            
            if mult_loc_check:
                mult_loc_check_overall = mult_loc_check
            
            ax.text(conseqs_x[1], y + conseqs_text_ydiff, prob_text, 
                    fontsize=text_font_size, fontweight=fweight, zorder=text_zorder, 
                    horizontalalignment='center', verticalalignment='center')
            
            if val_text:
                ax.text(conseqs_x[1], y - (conseqs_text_ydiff * 0.9), val_text, 
                        fontsize=text_font_size_2, fontweight=fweight, zorder=text_zorder, 
                        horizontalalignment='center', verticalalignment='top')
        
        if plot_legend:
            # The viewer will generally assume that the simulation starts at 
            # t = 0 years, but it is possible to start the simulation at another 
            # time.
            if time_array[0] != 0:
                time_label = 'Simulation from t$_0$ = {:.2f} years\nto t = {:.2f} years, '.format(
                    time_array[0], time_array[t_index])
                time_label += 'Results Shown\nfor t = {:.2f} years'.format(
                    time_array[t_index])
                prob_text = 'Probability of Exceeding Threshold\nBetween Times t$_0$ and t'
            else:
                time_label = 'Results Shown\nfor t = {:.2f} years'.format(
                    time_array[t_index])
                prob_text = 'Probability of Exceeding Threshold\nby t = {} years'.format(
                    time_array[t_index])
            
            ax.text(center_x - time_x_offset, center_y - time_y_offset, time_label, 
                    fontsize=text_font_size, fontweight=fweight, zorder=text_zorder, 
                    horizontalalignment='center', verticalalignment='center')
            
            analysis_type_text = ANALYSIS_TYPE_LABEL.get(analysis, None)
            
            analysis_value_adjustment = ''
            if not analysis == 'forward':
                analysis_value_adjustment = 'Mean '
                if realizations:
                    analysis_type_text = analysis_type_text.format(realizations)
                else:
                    analysis_type_text = None
            
            if analysis_type_text:
                ax.text(center_x - analysis_x_offset, center_y - analysis_y_offset, 
                        analysis_type_text, fontsize=text_font_size, 
                        fontweight=fweight, zorder=text_zorder, 
                        horizontalalignment='center', verticalalignment='center')
            
            ax.plot(legend_line_x, legend_line_y, color=whisker_color, 
                    linewidth=lwidth, linestyle=lstyle, zorder=whisker_zorder)
            
            ax.text(legend_line_x[1], legend_line_y[1] + prob_legend_y_diff, 
                    prob_text, fontsize=text_font_size_2, fontweight=fweight, 
                    zorder=text_zorder, horizontalalignment='center', 
                    verticalalignment='center')
            
            stand_dev_text = ''
            if analysis in ['lhs', 'parstudy']:
                stand_dev_text = '\u00B1 Stand. Dev. '
            
            value_text = '{}Value {}at Time = t'.format(
                analysis_value_adjustment, stand_dev_text)
            
            if mult_loc_check_overall:
                value_text += ',\nRanges Shown if Multiple Locations Are Used'
            
            ax.text(legend_line_x[1], legend_line_y[1] - val_legend_y_diff, 
                    value_text, fontsize=text_font_size_2, fontweight=fweight, 
                    zorder=text_zorder, horizontalalignment='center', verticalalignment='top')
            
            if plot_thresholds:
                ax.text(thresh_legend_x, thresh_legend_y, 
                        'Metric Threshold', fontsize=text_font_size, fontweight=fweight, 
                        color=thresh_text_color, zorder=text_zorder, 
                        horizontalalignment='center', verticalalignment='center')
                
                ax.text(expert_legend_x, expert_legend_y, 
                        'or Expert Review', fontsize=text_font_size, fontweight=fweight, 
                        color=expert_panel_color, zorder=text_zorder, 
                        horizontalalignment='center', verticalalignment='center')
        
        if '.' in name:
            name_main = name[0:name.index('.')]
            name_extension = name[name.index('.'):]
        else:
            name_main = name
            name_extension = '.png'
        
        if len(time_indices) > 1:
            name_main += '_tIndex_{:.0f}'.format(t_index)
        
        # Do not need the ticks
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        
        # Do not need the spline
        plt.setp(ax.spines.values(), alpha = 0)
        
        ax.set_aspect('equal')
        
        plt.savefig(os.sep.join([output_dir, name_main + name_extension]), 
                    dpi=figdpi)
        plt.close()


def get_result_text(prob, val, std_val, metric, show_std_option=True, 
                    analysis_type=None, unit=None):
    """
    Takes the probability, value, and standard deviation and produces the 
    formatted text used in make_bowtie_figure(). The variable mult_loc_check 
    indicates if multiple locations were used for the output type.
    """
    mult_loc_check = False
    
    if isinstance(prob, str):
        prob_text = prob
        val_text = None
    else:
        prob_text = '{:.2f} %'.format(prob * 100)
        
        # The value, std_vals, min_vals, and max_vals lists will have an entry 
        # for each location producing the output in question. Note that the min, 
        # max, and std values here reflect the differences across stochastic 
        # realizations - not across different locations. If the simulation is 
        # stochastic, the values being handled here in the val variable are the 
        # mean values across different realizations.
        if len(val) > 1:
            val_text = []
            
            min_val_diff_locs = np.min(val)
            max_val_diff_locs = np.max(val)
            
            if std_val:
                min_val_std = std_val[np.argmin(val)]
                max_val_std = std_val[np.argmax(val)]
            else:
                min_val_std = None
                max_val_std = None
            
            if (min_val_diff_locs == max_val_diff_locs) and (min_val_std == max_val_std):
                val_text = get_val_str(min_val_diff_locs, min_val_std, metric, 
                                       show_std_option=show_std_option, 
                                       analysis_type=analysis_type, unit=unit)
                
            else:
                min_val_str = get_val_str(min_val_diff_locs, min_val_std, metric, 
                                          show_std_option=show_std_option, 
                                          analysis_type=analysis_type, unit=unit)
                
                max_val_str = get_val_str(max_val_diff_locs, max_val_std, metric, 
                                          show_std_option=show_std_option, 
                                          analysis_type=analysis_type, unit=unit)
                
                val_text =  min_val_str + '\nto ' + max_val_str 
                
                mult_loc_check = True
        else:
            val_text = get_val_str(val[0], std_val[0], metric, 
                                   show_std_option=show_std_option, 
                                   analysis_type=analysis_type, unit=unit)
    
    return prob_text, val_text, mult_loc_check


def get_val_str(val, std_val, metric, show_std_option=True, analysis_type=None, 
                unit=None):
    """
    Takes a value, standard deviation, and metric type and produces a formatted 
    string reflecting the input. The metric type (e.g., brine_aquifer1 or 
    brine_aquifer2) is used to get the correct units from UNIT_DICT. Otherwise, 
    if a unit type is specified through the Contributors and Consequences .csv 
    files, that unit type is used.
    """
    if not analysis_type is None:
        if analysis_type == 'forward':
            show_std_option = False
        elif analysis_type in ['lhs', 'parstudy']:
            show_std_option = True
    
    if not unit is None:
        unit_type = unit
    elif metric in UNIT_DICT:
        unit_type = UNIT_DICT[metric]
    else:
        unit_type = ''
    
    # For the string to be formatted correctly, the unit type cannot be '$'
    if unit_type == '$':
        unit_type = 'dollars'
    
    if val == 0:
        val_text = '0'
    else:
        a, b = '{:.2e}'.format(val).split('e')
        b = int(b)
        
        val_text = r'${}\times10^{{{}}}$'.format(a, b)
    
    if std_val and show_std_option:
        val_text += ' \u00B1 '
        
        a, b = '{:.2e}'.format(std_val).split('e')
        b = int(b)
        
        if isinstance(metric, list):
            metric = metric[0]
        
        val_text += r'${}\times10^{{{}}}$'.format(a, b) + ' ' + unit_type
    else:
        val_text += ' ' + unit_type
    
    return val_text

