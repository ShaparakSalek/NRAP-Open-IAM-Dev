import os
import sys
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import openiam as iam
    import openiam.visualize as iam_vis
    from openiam import IAM_DIR
except ImportError as err:
    print('Unable to load openiam module: {}'.format(err))


STRAT_FIG_SIZE = (12, 10)
TIME_SERIES_FIG_SIZE = (13, 8)
AOR_FIG_SIZE = (7, 6)
TTFD_FIG_SIZE = (10, 8)
RADIAL_OBS_FIG_SIZE = (10, 8)

def process_plots(yaml_data, model_data, sm, s, output_list, analysis,
                  time_array, components, locations):
    """
    Analyze control file setup related to the plots setup to determine
    what kind of plots need to be produced.
    """

    output_dir = os.path.join(IAM_DIR, model_data['OutputDirectory'])

    plots = yaml_data['Plots']
    for p in plots:
        # Figure name
        save_filename = plots[p].get('savefig', p)  # plot name p is default file name
        save_filename = os.path.join(output_dir, save_filename)

        # Title of plot
        title = plots[p].get('Title', None)
        # Use or not subplots
        subplot = plots[p].get('subplot', {'use': False})
        if 'use' not in subplot:
            subplot['use'] = True

        if 'Data' in plots[p]:
            # Keep Data keyword for backward compatibility
            plots[p]['TimeSeries'] = plots[p]['Data']

        if 'TimeSeries' in plots[p]:
            iam_vis.time_series_plot(
                plots[p]['TimeSeries'], sm, s, plots[p], output_list,
                name=p, analysis=analysis, savefig=save_filename,
                title=title, subplot=subplot, plot_type=['real'],
                figsize=tuple(plots[p].get('figsize', TIME_SERIES_FIG_SIZE)))

        if 'TimeSeriesStats' in plots[p]:
            iam_vis.time_series_plot(
                plots[p]['TimeSeriesStats'], sm, s, plots[p], output_list,
                name=p, analysis=analysis, savefig=save_filename,
                title=title, subplot=subplot, plot_type=['stats'],
                figsize=tuple(plots[p].get('figsize', TIME_SERIES_FIG_SIZE)))

        if 'TimeSeriesAndStats' in plots[p]:
            iam_vis.time_series_plot(
                plots[p]['TimeSeriesAndStats'], sm, s, plots[p], output_list,
                name=p, analysis=analysis, savefig=save_filename,
                title=title, subplot=subplot, plot_type=['real', 'stats'],
                figsize=tuple(plots[p].get('figsize', TIME_SERIES_FIG_SIZE)))

        if 'AtmPlumeSingle' in plots[p]:
            satm = find_atm_comp(components)

            iam_vis.map_plume_plot_single(plots[p], p, sm, s, satm, time_array,
                                          output_dir, analysis=analysis)

        if 'AtmPlumeEnsemble' in plots[p]:
            satm = find_atm_comp(components)

            iam_vis.map_plume_plot_ensemble(plots[p], p, sm, s, satm, time_array,
                                            output_dir, analysis=analysis)

        if 'Stratigraphy' in plots[p]:
            iam_vis.stratigraphy_plot(
                yaml_data, model_data, sm, name=p,
                savefig=save_filename, title=title, figsize=tuple(plots[p].get(
                    'figsize', STRAT_FIG_SIZE)))

        if 'AoR' in plots[p]:
            iam_vis.area_of_review_plot(
                yaml_data, model_data, plots[p]['AoR'], sm, s, output_list,
                locations, name=p, analysis=analysis, savefig=save_filename,
                title=title, figsize=tuple(plots[p].get('figsize', AOR_FIG_SIZE)),
                save_results=True)

        if 'TTFD' in plots[p]:
            iam_vis.ttfd_plot(
                yaml_data, model_data, sm, s,
                output_dir, name=p, analysis=analysis,
                figsize=tuple(plots[p].get('figsize', TTFD_FIG_SIZE)),
                genfontsize=12, axislabelfontsize=14, titlefontsize=14,
                boldlabels=True)

        if 'GriddedRadialMetric' in plots[p]:
            iam_vis.gridded_radial_metric_plot(
                yaml_data, sm, output_dir, name=p,
                analysis=analysis, savefig=save_filename,
                figsize=tuple(plots[p].get('figsize', RADIAL_OBS_FIG_SIZE)),
                genfontsize=12, axislabelfontsize=14, titlefontsize=14,
                boldlabels=True)


def find_atm_comp(components):
    """ Return AtmosphericROM component if it was added to the system model."""

    for comp in components[::-1]:
        if isinstance(comp, iam.AtmosphericROM):
            atm_comp = comp
            break
    else:
        atm_comp = None
        logging.warning('Unable to find Atmospheric ROM for plume plots')

    return atm_comp
