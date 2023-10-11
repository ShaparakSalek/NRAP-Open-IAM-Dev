
# Default settings for TimeSeries plots
DEFAULT_TIME_SERIES_PLOT = 'TimeSeries'

DEFAULT_TIME_SERIES_SUBPLOT = {'Use': True}
DEFAULT_TIME_SERIES_USE_MARKERS = False
DEFAULT_TIME_SERIES_USE_LINES = True
DEFAULT_TIME_SERIES_VARY_LINES = False

LEAKAGEASSESSMENT_PLOT_NAMES = [{'Pressure_Plot': 'pressure',
                                 'CO2_Sat_Plot': 'CO2saturation',
                                 'CO2_Leakage_Plot': 'CO2_aquifer{}',
                                 'Brine_Leakage_Plot': 'brine_aquifer{}'}]


def set_up_leakage_assessment_workflow_plots(yaml_data, aquifer_name, figure_dpi=100,
                                             fig_size=None, extension=None):
    """
    Sets up the plot entry dictionary for the LeakageAssessment workflow.
    """
    plot_types_to_add = LEAKAGEASSESSMENT_PLOT_NAMES

    for plotNum, plot in enumerate(plot_types_to_add):

        for plotName in list(plot.keys()):

            if 'aquifer' in plot[plotName]:
                aquiferNum = aquifer_name[7:]

                metric = plot[plotName].format(aquiferNum)

            else:
                metric = plot[plotName]

            plot_type = yaml_data['Workflow']['Options'].get(
                'PlotType', DEFAULT_TIME_SERIES_PLOT)

            if 'subplot' in yaml_data['Workflow']['Options']:
                subplot = yaml_data['Workflow']['Options']['subplot']
            else:
                subplot = yaml_data['Workflow']['Options'].get(
                    'Subplot', DEFAULT_TIME_SERIES_SUBPLOT)

            use_markers = yaml_data['Workflow']['Options'].get(
                'UseMarkers', DEFAULT_TIME_SERIES_USE_MARKERS)

            use_lines = yaml_data['Workflow']['Options'].get(
                'UseLines', DEFAULT_TIME_SERIES_USE_LINES)

            vary_lines = yaml_data['Workflow']['Options'].get(
                'VaryLineStyles', DEFAULT_TIME_SERIES_VARY_LINES)

            yaml_data['Plots'][plotName + extension] = {
                plot_type: [metric], 'FigureDPI': figure_dpi,
                'Subplot': subplot, 'UseMarkers': use_markers,
                'UseLines': use_lines, 'VaryLineStyles': vary_lines}

            if fig_size:
                yaml_data['Plots'][plotName + extension]['FigureSize'] = fig_size

    return yaml_data
