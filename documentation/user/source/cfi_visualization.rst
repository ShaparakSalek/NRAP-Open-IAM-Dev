Setup of visualization options
------------------------------

The plot types available within NRAP-Open-IAM are: ``TimeSeries``, ``TimeSeriesStats``,
``TimeSeriesAndStats``, ``Stratigraphy``, ``AoR``, ``TTFD``, ``GriddedRadialMetric``,
``AtmPlumeSingle``, and ``AtmPlumeEnsemble``. We review the process of creating these
plots and then discuss each plot type separately.

To create a figure using a simulation run with a *.yaml* control file, the
the figure must be set up in the *Plots* section. Within the *Plots* section, the
user can have multiple entries that each represent a different plot to be created.
The names of these entries are defined by the user, and the names are generally used
as the file name for the corresponding figure. An exception occurs with some
plots. For example, with the ``TTFD`` plot type generates a large number of
plots, and the final file names will depend on the input used
(e.g., *TDS_Plume_Timings_Realization10.png*). The plot name provided can have
extensions appended that specify the file type of the resulting figure file
(e.g., *.png*, *.tiff*, or *.eps*). Below, we show two examples of ``TimeSeries`` plot
entries in a *.yaml* control file.

.. code-block:: python
   :lineno-start: 1

    Plots:
        Pressure_Figure:
            TimeSeries: [pressure]
        CO2_Sat_Figure.tiff:
            TimeSeries: [CO2saturation]

Since the first plot entry (*Pressure_Figure*) does not have an extension
(e.g., *.png* or *.tiff*) appended at the end (e.g., *Pressure_Figure.tiff*),
the produced figure will be of the default *.png* type. The second plot entry
(*CO2_Sat_Figure.tiff*), however, includes *.tiff* at the end of the name;
the resulting figure file will be of *.tiff* type. Note that if an extension
is provided when using a plot type that generates names (e.g., the ``TTFD`` plot type),
the generated names will still use the extension provided.

When we refer to an entry as being indented beneath another entry, we
mean that in a *.yaml* file the indented entry is on a lower line and preceeded
by additional four spaces. The indentation of one entry carries over to all entries
contained within that entry. In the example above, for example, *Pressure_Figure:*
is indented beneath ``Plots:``, and ``TimeSeries: [pressure]`` is indented beneath
*Pressure_Figure:*. ``TimeSeries`` is preceeded by eight spaces, while
*Pressure_Figure* is preceeded by four spaces. Additionally, the name of each entry
is followed by a colon (:). When entries have other entries indented beneath
them, the colon will be the last character on that line (e.g., ``Plots:`` or
``Pressure_Figure:``). When entries do not have other entries indented beneath
them and instead take an input value or list, then the colon is followed by that
input (e.g., ``TimeSeries: [pressure]``). Note that for the rest of this section,
we will not include the colon when discussing an entry; it is assumed
that each entry is followed by a colon.

All plot types have certain entries that are either required or
optional. Some of these optional entries are used by multiple plot types. To
avoid repeating the definitions of these entries, we first present the optional
entries used by multiple plot types. Then, we review each plot type.

Each plot type has the optional entry ``FigureDPI``.

* ``FigureDPI`` - the dots-per-inch (DPI) of the resulting figure(s) (default is
  100). Most figure types produce only one figure file, but plot types like
  ``Stratigraphy`` and ``TTFD`` can produce multiple figures from one entry.
  Larger DPIs will create high-resolution, high-quality figures, but the file
  sizes are also larger. File size and quality is also influenced by the extension
  used (e.g., *.png* or *.tiff*). Recommended ``FigureDPI`` values are between
  100 and 300.

Seven examples of ``FigureDPI`` entries in a *.yaml* control file are shown below.

.. code-block:: python
   :lineno-start: 1

    Plots:
        Stratigraphy_Figure.tiff:
            Stratigraphy:
                FigureDPI: 300
        Pressure_Figure:
            TimeSeriesStats: [pressure]
            FigureDPI: 100
        TDS_Volume_AoR.tiff:
            AoR: [Dissolved_salt_volume]
            FigureDPI: 200
        TTFD_Figures:
            TTFD:
                PlumeType: Dissolved_salt
                ComponentNameList: [GenericAquifer1, GenericAquifer2]
                FigureDPI: 200
        Gridded_Salt_Mass_Frac:
            GriddedRadialMetric:
                ComponentNameList: [GenericAquifer1]
                MetricName: Dissolved_salt_mass_fraction
                FigureDPI: 200
        Atmospheric_Plume_Single:
            AtmPlumeSingle:
                FigureDPI: 200
        Atmospheric_Plume_Probability.tiff:
            AtmPlumeEnsemble:
                FigureDPI: 300

Notice that the ``FigureDPI`` entries for the ``Stratigraphy``, ``TTFD``,
``GriddedRadialMetric``, ``AtmPlumeSingle``, and ``AtmPlumeEnsemble`` plots
are indented under the plot type. In contrast, the ``FigureDPI`` entries for the
``TimeSeriesStats`` and ``AoR`` plots are not indented beneath the plot type. This
discrepancy occurs because the ``TimeSeriesStats`` and ``AoR`` entries are followed
by a metric (e.g., [**pressure**]), while the other plot type entries are not.

The ``Stratigraphy``, ``TTFD``, ``GriddedRadialMetric``, ``AtmPlumeSingle``, and
``AtmPlumeEnsemble`` plot types all have the optional entries ``PlotInjectionSites``,
``InjectionCoordx``, ``InjectionCoordy``, ``SpecifyXandYLims``, and ``SaveCSVFiles``.

* ``PlotInjectionSites`` - an option to plot injection sites (default is ``False``).
  The only acceptable values are ``True`` or ``False``.

* ``InjectionCoordx`` - value or list of values for the x coordinate(s) of
  injection site(s) (default is None). The value(s) are in meters. This entry must
  be provided when using a ``LookupTableReservoir``, as that component type does
  not have an .injX attribute. Other reservoir types like ``SimpleReservoir`` or
  ``AnalyticalReservoir`` can be displayed without an InjectionCoordx entry.

* ``InjectionCoordy`` - value or list of values for the y coordinate(s) of
  injection site(s) (default is None). The value(s) are in meters. This entry must
  be provided when using a ``LookupTableReservoir``, as that component type does
  not have an .injY attribute. Other reservoir types like ``SimpleReservoir`` or
  ``AnalyticalReservoir`` can be displayed without an InjectionCoordy entry.

* ``SaveCSVFiles`` - an option to save results in *.csv* files. The only acceptable
  values are ``True`` or ``False``. The default value for ``AoR``, ``TTFD``,
  ``GriddedRadialMetric``, ``AtmPlumeSingle``, and ``AtmPlumeEnsemble`` plots is
  ``True``, while the default value for ``Stratigraphy`` plots is ``False``. For
  ``Stratigraphy`` plots, the *.csv* files contain unit thicknesses and depths
  across the domain. The *.csv* files are not saved by ``Stratigraphy`` plots when
  the simulation uses the ``LookupTableStratigraphy`` option.

If set up, ``SpecifyXandYLims`` is a dictionary containing two entries: ``xLims``
and ``yLims`` (i.e., ``xLims`` and ``yLims`` are indented beneath
``SpecifyXandYLims`` in a *.yaml* file).

* ``SpecifyXandYLims`` - a dictionary containing two optional entries related
  to the limits of the figure's x and y axes (default is None). Within
  this dictionary are the entries ``xLims`` and ``yLims``.

* ``xLims`` - an entry under ``SpecifyXandYLims`` containing a list of length two
  that represents the x-axis limits (e.g., ``xLims: [0, 1000]``; default is None).
  The values are in meters. The first and second values in the list are the
  lower and upper limits, respectively. If ``xLims`` is not provided or provided
  incorrectly, the figure will use the default approach for setting the
  x-axis limits.

* ``yLims`` - an entry under ``SpecifyXandYLims`` containing a list of length two
  that represents the y-axis limits (e.g., ``yLims: [0, 1000]``; default is None).
  The values are in meters. The first and second values in the list are the
  lower and upper limits, respectively. If ``yLims`` is not provided or provided
  incorrectly, the figure will use the default approach for setting the
  y-axis limits.

The ``Stratigraphy``, ``TTFD``, and ``AtmPlumeEnsemble`` plots also have the optional
entry ``SpecifyXandYGridLims``, which is a dictionary containing the ``gridXLims`` and
``gridYLims`` entries. ``AoR`` plots do not have grid entries because the x and y values
used are those of the ``OpenWellbore`` components. ``GriddedRadialMetric`` plots use the
radial grids produced by a component (e.g., a ``GenericAquifer`` component).

* ``SpecifyXandYGridLims`` - a dictionary containing two optional entries
  related to the x and y limits for the gridded data evaluated (default is None).
  In ``Stratigraphy`` plots, the gridded data are the three-dimensional planes
  depicting the the top of each unit. For ``TTFD`` and ``AtmPlumeEnsemble`` plots, the
  gridded data are the color-labelled values. Within this dictionary are the
  entries ``gridXLims`` and ``gridYLims``.

* ``gridXLims`` - an entry under ``SpecifyXandYGridLims`` containing a list of
  length two that represents the x-axis limits for the grid used to evaluate results
  (e.g., ``gridXLims: [100, 900]``; default is None). The values for ``gridXLims`` are
  in meters. The first and second values in the list are the lower and upper
  limits, respectively. If ``gridXLims`` is not provided or provided incorrectly,
  the figure will use the default approach for creating the gridded values.

* ``gridYLims`` - n entry under ``SpecifyXandYGridLims`` containing a list of
  length two that represents the y-axis limits for the grid used to evaluate results
  (e.g., ``gridYLims: [100, 900]``; default is None). The values for ``gridYLims`` are
  in meters. The first and second values in the list are the lower and upper
  limits, respectively. If ``gridYLims`` is not provided or provided incorrectly,
  the figure will use the default approach for creating the gridded values.

The ``Stratigraphy``, ``TTFD``, and ``AtmPlumeEnsemble`` plot types can all use
the optional entries ``xGridSpacing`` and ``yGridSpacing``:

* ``xGridSpacing`` - a horizontal distance (m) used as the interval between the
  grid points in the x-direction (default is None). If this entry is not setup,
  the x-coordinates of the grid points are defined using a default approach
  (1/100th of the range in x-values).

* ``yGridSpacing`` - a horizontal distance (m) used as the interval between the
  grid points in the y-direction (default is None). If this entry is not setup,
  the y-coordinates of the grid points are defined using a default approach
  (1/100th of the range in y-values).

The ``AoR`` and ``GriddedRadialMetric`` plot types both have the optional
entry ``TimeList``:

* ``TimeList`` - a list specifying the times (in years) for which to create separate
  figures (e.g., ``TimeList: [1, 5, 10]). Otherwise, one figure can be created for
  each timestep by having ``TimeList: All``. If TimeList is not entered for an ``AoR``
  plot, the figures created will show the maximum values for all locations across all
  model times. If ``TimeList`` is not entered for a ``GriddedRadialMetric`` plot,
  the default setting is ``TimeList: All``.

The ``TTFD`` and ``GriddedRadialMetric`` plot types both have the optional
entry ``ComponentNameList``:

* ``ComponentNameList`` - a list containing the names provided for each of the
  components producing output to be used for the creation of the figures (e.g.,
  ``ComponentNameList: [FutureGen2AZMI1, FutureGen2AZMI2]`` in
  *ControlFile_ex40.yaml*). Below, we show a section of the *.yaml* file for
  *ControlFile_ex40.yaml*. This section demonstrates where the name is provided
  for the *FutureGen2AZMI2* component. Below the excerpt is an example of how
  component names are set when using NRAP-Open-IAM in a script application.

Excerpt from *ControlFile_ex40* demonstrating how an aquifer component is given
the name FutureGen2AZMI2:

.. code-block:: python
   :lineno-start: 1

    FutureGen2AZMI2:
        Type: FutureGen2AZMI
        Connection: MultisegmentedWellbore1
        AquiferName: aquifer3
        Parameters:
            por: 0.132
            log_permh: -12.48
            log_aniso: 0.3
            rel_vol_frac_calcite: 0.1
        Outputs: [pH_volume, TDS_volume, Dissolved_CO2_volume,
                  Dissolved_CO2_dx, Dissolved_CO2_dy, Dissolved_CO2_dz]

Example of setting the component name (*FutureGen2AZMI2*) in a script application:

.. code-block:: python
   :lineno-start: 1

    fga = sm.add_component_model_object(FutureGen2AZMI(name='FutureGen2AZMI2', parent=sm))

The ``GriddedRadialMetric`` and ``AtmPlumeSingle`` plot types both have the optional
entry ``Realization``:

* ``Realization`` - the realization number for which to display results (default is 0).
  Note that this optional input is only used in simulations using Latin Hypercube Sampling
  (``lhs``) and Parameter Study (``parstudy``) analysis types. This input uses the indexing
  rules in Python, where 0 represents the first realization and (N - 1) represents the last
  (where N is the number of realizations).

Below are six plot entries in a *.yaml* file demonstrating the use of these
options. The ``PlumeType`` entry for the ``TTFD`` plot and the ``MetricName`` entry
for the ``GriddedRadialMetric`` plot are described further below.

.. code-block:: python
   :lineno-start: 1

    Plots:
        Strat_Plot.tiff:
            Stratigraphy:
                 FigureDPI: 300
                 PlotInjectionSites: True
                 PlotInjectionSiteLabels: True
                 InjectionCoordx: 200
                 InjectionCoordy: 200
                 SpecifyXandYLims:
                     xLims: [0, 400]
                     yLims: [0, 400]
                 SpecifyXandYGridLims:
                     gridXLims: [50, 350]
                     gridYLims: [50, 350]
                 xGridSpacing: 25
                 yGridSpacing: 25
        AoR_TDS:
            AoR: [Dissolved_salt_volume]
            FigureDPI: 200
            PlotInjectionSites: True
            InjectionCoordx: 200
            InjectionCoordy: 200
            SpecifyXandYLims:
                xLims: [0, 400]
                yLims: [0, 400]
            TimeList: All
        TTFD_TDS.tiff:
            TTFD:
                PlumeType: Dissolved_salt
                ComponentNameList: [GenericAquifer1, GenericAquifer2]
                FigureDPI: 300
                PlotInjectionSites: True
                InjectionCoordx: 200
                InjectionCoordy: 200
                SpecifyXandYLims:
                    xLims: [0, 400]
                    yLims: [0, 400]
                SpecifyXandYGridLims:
                    gridXLims: [50, 350]
                    gridYLims: [50, 350]
                 xGridSpacing: 5
                 yGridSpacing: 5
        Gridded_Salt_Mass_Frac.tiff:
            GriddedRadialMetric:
                MetricName: Dissolved_salt_mass_fraction
                ComponentNameList: [GenericAquifer1]
                FigureDPI: 300
                PlotInjectionSites: True
                InjectionCoordx: 200
                InjectionCoordy: 200
                TimeList: All
        ATM_single:
            AtmPlumeSingle:
                FigureDPI: 100
                PlotInjectionSites: True
                InjectionCoordx: 200
                InjectionCoordy: 200
                SpecifyXandYLims:
                    xLims: [-125, 325]
                    yLims: [-125, 325]
        ATM_Ensemble.tiff:
            AtmPlumeEnsemble:
                FigureDPI: 300
                PlotInjectionSites: True
                InjectionCoordx: 200
                InjectionCoordy: 200
                PlotReceptors: False
                SpecifyXandYLims:
                    xLims: [-125, 325]
                    yLims: [-125, 325]
                SpecifyXandYGridLims:
                    gridXLims: [-100, 300]
                    gridYLims: [-100, 300]
                xGridSpacing: 1
                yGridSpacing: 1

TimeSeries, TimeSeriesStats, and TimeSeriesAndStats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``TimeSeries``, ``TimeSeriesStats``, and ``TimeSeriesAndStats`` plot types
are used to display results varying over time. Although this section
covers three plot types, these plot types are different variations of
the same type of plot.

``TimeSeries`` plots are line plots of results varying over time. The number
of lines in the resulting figure depends on the setup of the scenario. For example,
components and associated locations entered in the *.yaml* file can define the
number of curves shown in the figure but only the components that produce the metric
being plotted (e.g., **pressure** or **brine_aquifer1**) influence the number
of lines created for that particular metric.

``TimeSeriesStats`` and ``TimeSeriesAndStats`` plots can only be produced for simulations
using the Latin Hypercube Sampling (LHS, ``lhs`` in a control file setup)
or Parameter Study (``parstudy`` in a control file setup) analysis
types (not the ``forward`` analysis type). Simulations using the ``lhs`` and ``parstudy``
analysis types create separate simulations (i.e., different realizations) that
explore the parameter space. The parameters varied are those entered with
minimum and maximum values, which are meant to model uniform distribution. Consider,
for example, a case of ``TimeSeriesStats`` plot for an LHS run with 30
realizations. The ``ModelParams`` section of the *.yaml* file would be similar to
this excerpt from *ControlFile_ex4a.yaml*:

.. code-block:: python
   :lineno-start: 1

    ModelParams:
        EndTime: 10
        TimeStep: 1.0
        Analysis:
            Type: lhs
            siz: 30
        Components: [SimpleReservoir1,
                     OpenWellbore1,
                     CarbonateAquifer1]
        OutputDirectory: output/output_ex4a_{datetime}
        Logging: Debug

The entries ``Type: lhs`` and ``siz: 30`` under ``Analysis`` specify the run as an
LHS simulation with 30 realizations. Each realization will use different values
for the parameters that are setup to vary. In a ``TimeSeries`` plot, for each realization
a corresponding metric (e.g., **pressure**) will be represented by a separate line.

If an LHS or parstudy simulation uses many realizations and many component
locations, the ``TimeSeries`` plot could become visually unclear. To avoid a lack of visual
clarity, ``TimeSeriesStats`` plots show the basic information about the distribution
of results. The plot produces lines representing mean and median values as well
as shaded regions showing the four quartiles of the distribution varying over time
(0th to 25th, 25th to 50th, 50th to 75th and 75th to 100th percentiles).

``TimeSeriesAndStats`` plots combine the approaches of ``TimeSeries``
and ``TimeSeriesStats`` plots. The mean, median, and quartiles are shown along
with line graphs for each realization.

``TimeSeries`` and ``TimeSeriesAndStats`` plots can have several optional
entries: ``UseMarkers``, ``VaryLineStyles``, ``UseLines``, and ``FigureDPI``
(``FigureDPI`` is described above).

* ``UseMarkers`` - an option to show results with values annotated with markers
  like circles and squares (default is ``False``). The only acceptable values
  are ``True`` or ``False``. If markers are used, the colors of markers and lines
  will vary in the normal manner (i.e., a rotation through the default
  matplotlib color order).

* ``VaryLineStyles`` - an option to vary the line styles used (default is ``False``).
  The only acceptable values are ``True`` or ``False``. The matplotlib line styles
  used are 'solid', 'dotted', 'dashed', and 'dashdot'. Line colors will still
  vary in the normal manner.

* ``UseLines`` - an option to show results with lines (default is ``True``). The only
  acceptable values are ``True`` or ``False``. If neither markers nor lines are used,
  the plot will not show any results. One should only set ``UseLines`` to ``False``
  if ``UseMarkers`` is set to ``True``. If ``UseLines`` is set to ``False``,
  ``VaryLineStyles``  will automatically be set to ``False``, regardless
  of the entry provided in the *.yaml* file.

These entries are not indented under ``TimeSeries`` or ``TimeSeriesAndStats`` in a
*.yaml* file, but are instead indented under the figure name. For example, these
entries are used in *ControlFile_ex1b.yaml*:

.. code-block:: python
   :lineno-start: 1

    Plots:
        Pressure_plot_v1:
            TimeSeries: [pressure]
            UseMarkers: False
            UseLines: True
            VaryLineStyles: True
            FigureDPI: 150
        Pressure_plot_v2.tiff:
            TimeSeries: [pressure]
            UseMarkers: True
            UseLines: False
            VaryLineStyles: False
            FigureDPI: 400

If ``UseMarkers``, ``VaryLineStyles``, and ``UseLines`` are provided for a
``TimeSeriesStats`` plot, the entries will have no effect (i.e., they do not
influence the mean and median lines or the shaded quartiles).

For examples of ``TimeSeries`` plots, see control file examples 1a, 1b, 2, 3, 7a,
7b, and 14. For examples of ``TimeSeriesStats`` plots, see control file examples
4a, 4b, 6, 8, 15, and 39. For examples of ``TimeSeriesAndStats`` plots, see control
file examples 4a, 14, and 40.

Stratigraphy
~~~~~~~~~~~~
``Stratigraphy`` plots are three-dimensional figures showing the specified
stratigraphy as well as features like wellbores and injection sites. These plots
can vary with the approach used for the stratigraphy. For example, a ``strike`` and
``dip`` can be assigned in the ``Stratigraphy`` section of a *.yaml* control file.
Alternatively, the ``LookupTableStratigraphy`` option allows one to create the
domain's stratigraphy with a *.csv* file containing unit thicknesses. ``Stratigraphy``
plots also work for simulations with spatially uniform unit thicknesses.

First, we discuss the use of a ``strike`` and ``dip`` options. The ``Stratigraphy``
section from *ControlFile_ex33.yaml* is shown below:

.. code-block:: python
   :lineno-start: 1

    Stratigraphy:
        spatiallyVariable:
            strikeAndDip:
                strike: 315
                dip: 5
                dipDirection: NE
                coordxRefPoint: 1200
                coordyRefPoint: 1200
        numberOfShaleLayers:
            vary: False
            value: 3
        shale1Thickness:
            value: 750.0
            vary: False
        shale2Thickness:
            value: 950.0
            vary: False
        shale3Thickness:
            value: 200
            vary: False
        aquifer1Thickness:
            vary: False
            value: 200
        aquifer2Thickness:
            vary: False
            value: 200
        reservoirThickness:
            vary: False
            value: 150

To set up spatially variable stratigraphy, one can use ``spatiallyVariable`` keyword
indented under ``Stratigraphy``. To use strike and dip values, the ``strikeAndDip`` keyword
needs to be indented under ``spatiallyVariable``. The entries indented under ``strikeAndDip``
are as follows:

* ``strike`` - the strike of the units in degrees clockwise from north in a map
  view presentation. For example, strike values of 0 or 180 make the units
  strike north/south; strike values of 90 or 270 make the units strike
  east/west, and strike values of 30 or 210 make the units strike
  northeast/southwest. Acceptable values are in a range between 0 to 360.

* ``dip`` - the dip of the units in degrees, where a positive value corresponds
  with unit depths increasing in the ``dipDirection`` provided. Acceptable values
  range from 0 to less than 90.

* ``dipDirection`` - the dip direction provided in a cardinal direction -
  N, E, S, W, NE, SE, SW, or NW. Note that this entry must be compatible with
  the ``strike`` entry. For example, units cannot strike north/south and dip to
  the north, but they could strike north/south and dip to the east or west.

* ``coordxRefPoint`` - the x-coordinate (m) of the reference point. The unit
  thicknesses provided for the reference point are used to calculate unit
  thicknesses across the domain.

* ``coordyRefPoint`` - the y-coordinate (m) of the reference point. The unit
  thicknesses provided for the reference point are used to calculate unit
  thicknesses across the domain.

Note that the unit thicknesses indented under ``Stratigraphy`` are those at the
reference point (x = ``coordxRefPoint``, y = ``coordyRefPoint``). When using the
``strikeAndDip`` option, unit thicknesses in other parts of the domain are
calculated in relation to this reference point. Other ``Stratigraphy`` component
parameters like *numberOfShaleLayers* and *datumPressure* cannot vary across the
domain. Note that units can effectively pinch out, although the thicknesses will
only be reduced to the minimum value of 1 m. Additionally, while the ``strike``
and ``dip`` option will make some units thicker (e.g., increasing the thickness
of the the top shale so that the units beneath it have greater depths), each unit
thickness cannot exceed the maximum value of 1600 m.

To use the ``LookupTableStratigraphy`` approach, one can use ``spatiallyVariable`` indented
under ``Stratigraphy`` and then ``LookupTableStratigraphy`` keyword indented under
``spatiallyVariable``. This approach is demonstrated in *ControlFile_ex38.yaml*:

.. code-block:: python
   :lineno-start: 1

    Stratigraphy:
        spatiallyVariable:
            LookupTableStratigraphy:
                FileName: 'stratigraphy.csv'
                FileDirectory: 'examples/Control_Files/input_data/ex38'
                MaxPointDistance: 100

The entries indented under ``LookupTableStratigraphy`` are as follows:

* ``FileName`` - the name of the *.csv* file containing unit thicknesses and other
  ``Stratigraphy`` component parameters (*numberOfShaleLayers*, *datumPressure*, and
  *depth*).

* ``FileDirectory`` - the directory containing the *.csv* file referenced by
  *FileName*. The directory is given in relation to the main directory used for
  the NRAP-Open-IAM installation being used but ``FileDirectory`` can also provide
  an entire path name like

    C:\Users\UserName\Documents\NRAPOpenIAM\examples\Control_Files\input_data\ex38.

* ``MaxPointDistance`` - to set unit thicknesses at each location evaluated in
  the domain, each location must be within a certain distance of a point in
  the *.csv* file referenced with ``FileName``. ``MaxPointDistance`` is that maximum
  distance (m) (default is 100 m). If a location in the domain is not close
  enough to a point in the *.csv* file, the simulation will return an error.
  Users can avoid this error by setting ``MaxPointDistance`` to a higher value,
  but using too high a value could lead to inaccurate depictions of the
  domain's stratigraphy. ``MaxPointDistance`` is intended to help ensure that
  ``LookupTableStratigraphy`` *.csv* files include sufficient information. It is the
  user's responsibility to make sure that the *.csv* file contains sufficient
  information and the ``MaxPointDistance`` is not too high.

The first two columns of a ``LookupTableStratigraphy`` *.csv* file are x and y
coordinates (m) with the columns named 'x' and 'y', respectively.
Any unit thicknesses (m) that vary with x and y values should be listed in
columns with the same number of rows as the x and y columns. The thicknesses
specified in a particular row of the *.csv* file correspond to the x and y values
from the same row. If a unit thickness does not vary with x and y values,
that unit thickness can be displayed in a column with a single row. A location in
the domain will be assigned the unit thicknesses from the closest location in the
``LookupTableStratigraphy`` *.csv* file - if the closest location is within
``MaxPointDistance`` of the location. For an example, see the *stratigraphy.csv*
file in the directory *examples/Control_Files/input_data/ex38*.

Note that ``Stratigraphy`` plots created for simulations using ``LookupTableStratigraphy``
will not have three-dimensional planes. Instead, the tops of each unit are plotted
as squares along each wellbore.

``Stratigraphy`` plots can have several optional entries: ``PlotWellbores``,
``PlotWellLabels``, ``PlotInjectionSites``, ``PlotInjectionSiteLabels``,
``InjectionCoordx``, ``InjectionCoordy``, ``PlotStratComponents``,
``StrikeAndDipSymbol``, ``SpecifyXandYLims``, ``SpecifyXandYGridLims``,
``xGridSpacing``, ``yGridSpacing``, ``View``, ``SaveCSVFiles``, and ``FigureDPI``.
Four of these entries (``StrikeAndDipSymbol``, ``SpecifyXandYLims``,
``SpecifyXandYGridLims``, and ``View``) are dictionaries containing additional
entries (i.e., more entries indented beneath them in a *.yaml* file). The entries
``SpecifyXandYLims``, ``SpecifyXandYGridLims``, ``xGridSpacing``, ``yGridSpacing``,
``SaveCSVFiles``, ``PlotInjectionSites``, ``InjectionCoordx``, and
``InjectionCoordy`` are described above.

* ``PlotWellbores`` - an option to plot wellbores as vertical lines (default is
  ``True``). The only acceptable values are ``True`` or ``False``.

* ``PlotWellLabels`` - an option to show text labels specifying wellbore types
  and numbers (default is ``True``). For example, the labels could be "Open Wellbore
  1" for an Open Wellbore, "M.S. Wellbore 1" for a MultiSegmented Wellbore,
  or "Cemented Wellbore 1" for a Cemented Wellbore. The only acceptable values
  are ``True`` or ``False``.

* ``PlotInjectionSiteLabels`` - an option to show a text label for the injection
  site(s) (default is ``False``).

* ``PlotStratComponents`` - the option to plot squares along each wellbore at
  the depths at which the wellbore intersects the top of a unit (default is ``False``).
  The tops of shales are shown with red squares, while the tops of aquifers
  are shown with blue squares. The only acceptable values are ``True`` or ``False``.

* ``StrikeAndDipSymbol`` - a dictionary containing four optional entries related
  to the strike and dip symbol shown in the figure (default is None). Within
  this dictionary are the entries ``PlotSymbol``, ``coordx, coordy, and length.

* ``PlotSymbol`` - an entry under ``StrikeAndDipSymbol`` that specifies whether to
  show the strike and dip symbol (default is ``True``). The only acceptable values
  are ``True`` or ``False``.

* ``coordx`` - an entry under ``StrikeAndDipSymbol`` that specifies the x-coordinate
  at which to plot the strike and dip symbol (default is None). If ``coordx`` is
  not setup, the graph will use a default location (which depends on the domain).

* ``coordy`` - an entry under ``StrikeAndDipSymbol`` that specifies the y-coordinate
  at which to plot the strike and dip symbol (default is None). If ``coordy`` is
  not setup, the graph will use a default location (which depends on the domain).

* ``length`` - an entry under ``StrikeAndDipSymbol`` that specifies the length scale
  (m) of the strike and dip symbol (default is None). For flat-lying units, the
  length is the diameter of the circular symbol used. For dipping units, the
  length applies to the line going in direction of strike (not the line in
  the dip direction). If ``length`` is not provided, the graph will use a
  calculated length (which depends on the domain).

* ``View`` - a dictionary containing two optional entries related to the
  perspective of the three-dimensional graph (default is None). Within this
  dictionary are the entries ``ViewAngleElevation`` and ``ViewAngleAzimuth``.
  A separate version of the figure is created for each combination of
  the ``ViewAngleElevation`` and ``ViewAngleElevation`` entries, where
  the first values in the keywords list are used for the same graph and so on.

* ``ViewAngleElevation`` - an entry under ``View`` containing a list of the
  elevation angles (in degrees) to use in the ``Stratigraphy`` plot(s) (default is
  [10, 30]). Values must be between -90 and 90. See the matplotlib
  documentation regarding view angles. This list must have the same length as
  the ``ViewAngleAzimuth`` list.

* ``ViewAngleAzimuth`` - an entry under ``View`` containing a list of the
  azimuth angles (in degrees) to use in the ``Stratigraphy`` plot(s) (default is
  [10, 30]). Values must be between 0 and 360. See the matplotlib
  documentation regarding view angles. This list must have the same length as
  the ``ViewAngleElevation`` list.

Two examples of *.yaml* entries for ``Stratigraphy`` plots are shown below. The first
entry uses the default settings, while the second entry specifies each option.
Since the simulation uses a ``LookupTableReservoir``, the entry has to include
``InjectionCoordx`` and ``InjectionCoordy``. ``InjectionCoordx`` and
``InjectionCoordy`` are not required when using another type of reservoir
component with option ``PlotInjectionSites: True``.

.. code-block:: python
   :lineno-start: 1

    Plots:
        Strat_Plot_Default_Settings:
            Stratigraphy:
        Strat_Plot.tiff:
            Stratigraphy:
                 FigureDPI: 500
                 PlotInjectionSites: True
                 PlotInjectionSiteLabels: True
                 InjectionCoordx: 200
                 InjectionCoordy: 200
                 PlotWellbores: True
                 PlotWellLabels: True
                 PlotStratComponents: True
                 SaveCSVFiles: False
                 SpecifyXandYLims:
                     xLims: [0, 400]
                     yLims: [0, 400]
                 SpecifyXandYGridLims:
                     gridXLims: [25, 375]
                     gridYLims: [25, 375]
                 StrikeAndDipSymbol:
                     PlotSymbol: True
                     coordx: 100
                     coordy: 300
                     length: 75
                 View:
                     ViewAngleElevation: [5, 10, 5, 10]
                     ViewAngleAzimuth: [300, 300, 310, 310]

For examples of ``Stratigraphy`` plots, see examples *ControlFile_ex33.yaml*-*ControlFile_ex38.yaml*.
For examples of using ``Stratigraphy`` plots in a script application, see files
*iam_sys_reservoir_mswell_stratplot_dipping_strata.py* and
*iam_sys_reservoir_mswell_stratplot_no_dip.py*.

AoR
~~~
Area of Review (``AoR``) plots are developed to estimate the AoR needed for a geologic
carbon storage project based on the spatial extent of reservoir impacts (pressure
and |CO2| saturation) and potential aquifer impacts (dissolved salt and dissolved
|CO2| plume volumes). The potential extent is found by distributing ``OpenWellbore``
components across the domain. We recommend setting ``OpenWellbore`` locations
using the grid placement option (see examples *ControlFile_ex31a.yaml*  to *ControlFile_ex31d.yaml*).
The ``OpenWellbore`` (components) are hypothetical and used to consider the aquifer impacts
that could occur if a leakage pathway (extending from the reservoir to the aquifer being
considered) was available at each ``OpenWellbore`` location. The approach used for ``AoR``
plots is based on the work :cite:`BACON2020`.

Note that the ``AoR`` plot type is meant to be used only for one aquifer at a time,
with that aquifer being represented by only one type of aquifer component
(e.g., representing contaminant spread in aquifer 2 with a ``FutureGen2Aquifer``
component). For example, file *ControlFile_ex31a.yaml* has ``SimpleReservoir``
components that provide the input for ``OpenWellbore`` components, and the ``OpenWellbore``
components provide input to ``FutureGen2Aquifer`` components. The ``FutureGen2Aquifer``
components are set up to represent aquifer 2. If the user added an entry to the *.yaml*
file for a ``FutureGen2AZMI`` aquifer component representing aquifer 1, the ``AoR`` plot
could not make plots representing the impacts on both aquifers 1 and 2. In this
case, one would need to create a separate *.yaml* file that creates ``AoR`` plots just
for aquifer 1.

``AoR`` plots can be created for the following types of outputs: **pressure**,
**CO2saturation**, **pH_volume**, **TDS_volume**, **Dissolved_CO2_volume**,
and **Dissolved_salt_volume**. The ``AoR`` plot type examines these metrics
at each location in the domain (i.e., each hypothetical ``OpenWellbore``
location) and displays the maximum value over time (across all times or at specific
times, depending on the ``TimeList`` entry provided; this entry is discussed below).
For LHS simulations, the ``AoR`` plot displays the maximum values over time at
each location from all LHS realizations. This approach is meant to depict how
severe the reservoir and aquifer impacts could become. Using the ``AoR`` plot type
leads to the creation of *.csv* files containing the values shown in the ``AoR`` plots.
Note that model run times can increase dramatically with the number of ``OpenWellbore``
locations. Additionally, some aquifer components generally require longer model run
times (e.g., ``GenericAquifer``) in comparison with other aquifer components
(e.g., ``FutureGen2Aquifer``). Also note that ``FutureGen2Aquifer`` is meant to be
set up for aquifers with bottom depths <= 700 m, while ``FutureGen2AZMI`` is meant
to be set up for aquifers with bottom depths >= 700 m.

When using the ``AoR`` plot type, we recommend setting ``GenerateOutputFiles`` and
``GenerateCombOutputFile`` to ``False`` in the ``ModelParams`` section of the *.yaml* file.
The large number of ``OpenWellbore`` locations commonly used for ``AoR`` plots causes
a large number of output files. A reservoir and aquifer component is created for
each ``OpenWellbore`` location, and every component will have its output saved. The
``.csv`` files created for the ``AoR`` plots contain all of the necessary information
and these files are much smaller in size.

``AoR`` plots can have six optional entries: ``PlotInjectionSites``, ``InjectionCoordx``,
``InjectionCoordy``, ``SaveCSVFiles``, ``FigureDPI``, and ``TimeList``. All of these
entries are described above.

If the ``TimeList`` entry is not provided for an ``AoR`` plot, the figure will show the
maximum values at each location across all model times. If ``TimeList`` is provided
as a list of times in years (e.g., ``TimeList: [1, 5, 10]`` or ``TimeList: [10]``),
then the figures created will represent the maximum values at each location at the
specified time(s). Otherwise, an AoR figure can be made for every model time by providing
``TimeList: All``. Evaluating how the potential impacts of a project change over time
can inform, for example, how the required extents of surveying efforts change
over time (i.e., discovering and effectively plugging legacy wells at larger distances
from the injection site).

Below is an example of two ``AoR`` plot entries in a *.yaml* file. The first entry
uses the default settings, while the second specifies all available options.
Since the simulation uses a ``LookupTableReservoir`` this example includes
``InjectionCoordx`` and ``InjectionCoordy``. These inputs are not required
for other reservoir component types.

.. code-block:: python
   :lineno-start: 1

    Plots:
        AoR_pH_Default_Settings:
            AoR: [pH_volume]
        AoR_TDS.tiff:
            AoR: [TDS_volume]
            PlotInjectionSites: True
            InjectionCoordx: 2.37e5
            InjectionCoordy: 4.41e6
            FigureDPI: 300
            SaveCSVFiles: False
            TimeList: All

For examples of AoR plots, see *ControlFile_ex31a.yaml* to *ControlFile_ex32c.yaml*.

TTFD
~~~~
Using the time to first detection (``TTFD``) plot type can produce three types of
figures: maps of earliest plume timings across the domain (i.e., the earliest
time at which the plume type occurs in each part of the aquifer(s) considered),
maps showing the ``TTFD`` provided by the entered monitoring locations, and maps of
the probability of plume occurrence in the aquifer(s) considered. The figures
with the ``TTFD`` from monitoring locations are only created if monitoring locations
are entered. The maps of plume probabilities are only created if the analysis type
is Latin Hypercube Sampling (``lhs``) or Parameter Study (``parstudy``). Note that plume
probabilities are calculated as the number of realizations in which a plume
occurred at each location divided by the total number of realizations.

The ``TTFD`` plot type requires the use of at least one of the following aquifer
component types (with the component(s) set up to represent the aquifer(s)
considered): ``CarbonateAquifer``, ``FutureGen2Aquifer``, ``FutureGen2AZMI``, ``GenericAquifer``,
``DeepAlluviumAquifer``, or ``DeepAlluviumAquiferML``. Note that the ``FutureGen2Aquifer``
component is used for aquifers with bottom depths <= 700 m, while the ``FutureGen2AZMI``
component is used for aquifers with bottom depths >= 700 m. The aquifer component(s)
must also produce the plume dimension metrics associated with the plume type
considered (e.g., **TDS_dx**, **TDS_dy**, and **TDS_dz** for TDS plumes). Note that
``CarbonateAquifer`` components do not produce plume dimension outputs for different
plume types, so the required outputs when using ``CarbonateAquifer`` are **dx** and **dy**
(which represent the lengths of the impacted aquifer volume in the x-
and y-directions, respectively).

The plume timing and plume probability figures made with the ``TTFD`` plot type show
four subplots. Each subplot contains a quarter of the depth range from the
top of the reservoir to the surface. Each subplot contains the results for
sections of aquifers within the corresponding depth range. If monitoring sensor
locations are provided, each subplot will also show any sensors with depth (z) values
in the subplot's depth range as black triangles. Because there are multiple z grid
points within each subplot, there can be different layers of results displayed.
The code is set up to make the top layer shown be the layer with the lowest
plume timing or highest plume probability (for the corresponding figure types).
The matplotlib function used to display results by color (contourf) can fail to
display results when there are very few points with results in a layer. To
address such situations, if there are fewer than 25 points with results we
display each value as a color-labelled circle.

While the plume timing plots show the earliest plume timings at each grid location
across the domain, the monitoring ``TTFD`` plots only display plume timings that are
sufficiently close to the sensor location(s) provided. The purpose of such graphs
is to show when the sensors used could warn site operators that an aquifer has
been impacted. If the chosen sensor x, y, and z values do not provide any warning
of plumes in an aquifer, and there are plumes in that aquifer, then the monitoring
locations should be changed. The distance over which sensors can detect a plume
are controlled by the ``VerticalWindow`` and ``HorizontalWindow`` entries, which are
discussed below. Note that the ``TTFD`` plot type can produce output for the DREAM
tool (Design for Risk Evaluation And Management) if ``WriteDreamOutput`` is set to
``True`` (see below). DREAM is designed to optimize the placement of monitoring
sensors.

Unlike most other plot types, the ``TTFD`` plot type has two required entries:
``PlumeType`` and ``ComponentNameList``. ``TTFD`` plots will not be produced
without appropriate input for these entries. ``ComponentNameList`` is discussed
above.

* ``PlumeType`` - the type of plume metric being considered. Acceptable values
  are *Pressure*, *pH*, *TDS*, *Dissolved_CO2*, *Dissolved_salt*, and *CarbonateAquifer*.
  The dx, dy, and dz metrics (e.g., **Dissolved_CO2_dz**) for the PlumeType used
  must be produced by the aquifer components listed in ``ComponentNameList``. The
  dz metrics are not required when using ``CarbonateAquifer`` components, however,
  as these components do not produce a dz plume metric. Additionally, when
  using ``PlumeType: CarbonateAquifer`` the plume timing and plume probability
  figures do not have different subplots for different depth ranges.

The ``TTFD`` plot type can have the following optional entries: ``MonitoringLocations``,
``SaveCSVFiles``, ``WriteDreamOutput``, ``SpecifyXandYLims``, ``NumZPointsWithinAquifers``,
``NumZPointsWithinShales``, ``xGridSpacing``, ``yGridSpacing``, ``SpecifyXandYGridLims``,
``PlotInjectionSites``, ``InjectionCoordx``, ``InjectionCoordy``, and ``FigureDPI``.
Three of these entries (``MonitoringLocations``, ``SpecifyXandYLims``, and
``SpecifyXandYGridLims``) are dictionaries containing additional entries
(i.e., entries indented beneath mentioned keywords in a *.yaml* file).
All of these entries except for ``MonitoringLocations``, ``WriteDreamOutput``,
``NumZPointsWithinAquifers``, and ``NumZPointsWithinShales`` are described above.

The ``NumZPointsWithinAquifers``, ``NumZPointsWithinShales``, ``xGridSpacing``,
``yGridSpacing``, and ``SpecifyXandYGridLims`` entries all relate to the x-, y-,
and z-coordinates of the grids used to evaluate plume extents and timings.
The dx, dy, and dz plume dimension metrics (e.g., *pH_dy* or *TDS_dz*) are used
to evaluate whether each (x, y, z) of a grid is within a plume area for
each model timestep. Note that ``NumZPointsWithinAquifers`` and
``NumZPointsWithinShales`` do not have an effect when setup
``PlumeType: CarbonateAquifer`` is used because that ``CarbonateAquifer``
component does not produce a dz plume metric.

* ``MonitoringLocations`` - a dictionary containing five optional entries related
  to the sensors used to detect aquifer impacts. The five optional entries are
  ``coordx``, ``coordy``, ``coordz``, ``HorizontalWindow``, and ``VerticalWindow``.
  Note that the lists provided for ``coordx``, ``coordy``, and ``coordz`` must all
  have the same length (although ``coordz`` is not used with option
  ``PlumeType: CarbonateAquifer``).

* ``coordx`` - an entry under ``MonitoringLocations`` that specifies the
  x-coordinate(s) (m) of monitoring sensor(s), if any sensors are used. This entry
  must be provided as a list, even if only one location is used (e.g., [100]
  or [100, 200]).

* ``coordy`` - an entry under ``MonitoringLocations`` that specifies the
  y-coordinate(s) (m) of monitoring sensor(s), if any sensors are used. This entry
  must be provided as a list, even if only one location is used (e.g., [100]
  or [100, 200]).

* ``coordz`` - an entry under ``MonitoringLocations`` that specifies the depth(s)
  (z-coordinate(s), (m)) of monitoring sensor(s), if any sensors are used. Note that
  for this entry, depths beneath the surface are taken as negative values.
  This entry must be provided as a list, even if only one location is used
  (e.g., [-500] or [-500, -400]). The ``coordz`` entry is not required when using
  an option ``plumeType: CarbonateAquifer``, as the ``CarbonateAquifer``
  component does not produce a dz plume metric.

* ``HorizontalWindow`` - a (maximum) horizontal distance (m) from which monitoring
  sensor(s) will detect plumes (default is 1). For example, if the HorizontalWindow
  is 5 m, then the sensor will detect any plume at grid locations within 5 m
  of the sensor's ``coordx`` and ``coordy`` values (if the plume is also within
  ``VerticalWindow`` of the sensor's ``coordz`` value). This entry is meant to represent
  the sensitivity of a sensor, but that consideration must also involve the
  threshold used for the plume type considered (if the aquifer component has
  a user-defined threshold for plume detection). For example, **Dissolved_salt**
  plumes from the ``GenericAquifer`` are influenced by the **dissolved_salt_threshold**
  parameter. In contrast, the ``FutureGen2Aquifer`` component defines TDS plumes
  where the relative change in TDS is > 10% (i.e., no user-defined threshold).
  The inclusion of plumes at nearby grid points is also dependent on the spacing
  of grid points; the x- and y-spacings are controlled by ``xGridSpacing`` and
  ``yGridSpacing``, while the z-spacing is controlled by ``NumZPointsWithinAquifers``
  and ``NumZPointsWithinShales``. Note that the grid is made to include the x-, y-,
  and z-coordinates for monitoring locations, so there will always be a grid point
  for each monitoring sensor.

* ``VerticalWindow`` - a (maximum) vertical distance (m) from which monitoring
  sensor(s) will detect plumes (default is 1). For example, if the ``VerticalWindow``
  is 5 m, then the sensor will detect any plume within 5 m of the sensor's
  ``coordz`` values (if the plume is also within ``HorizontalWindow`` of the
  sensor's ``coordx`` and ``coordy`` value). This entry is meant to represent the
  sensitivity of a sensor, but that consideration must also involve the threshold
  used for the plume type considered (if the aquifer component has a user-defined
  threshold for plume detection). For example, **Dissolved_CO2** plumes from the
  ``GenericAquifer`` are influenced by the **dissolved_co2_threshold** parameter. In
  contrast, the ``FutureGen2Aquifer`` component defines pH plumes where the
  absolute change in pH is > 0.2 (i.e., no user-defined threshold). The
  inclusion of plumes at nearby grid points is dependent on the spacing of
  grid points; the x- and y-spacings are controlled by ``xGridSpacing`` and
  ``yGridSpacing``, while the z-spacing is controlled by ``NumZPointsWithinAquifers``
  and ``NumZPointsWithinShales``. Note that the grid is made to include the x-, y-,
  and z-coordinates for monitoring locations, so there will always be a grid point
  for each monitoring sensor.

* ``WriteDreamOutput`` - the option to create *.iam* files containing plume timing
  results (default is ``False``). These *.iam* files are the input for the DREAM
  program. DREAM is the Design for Risk Evaluation And Management tool, which
  was also developed by NRAP. The only acceptable values are ``True`` or ``False``.

* ``NumZPointsWithinAquifers`` - the number of z-grid points extending from the
  bottom to the top of each aquifer (default is 10). The points are equally
  spaced.

* ``NumZPointsWithinShales`` - the number of z-grid points extending from the
  bottom to the top of each shale (default is 3). The points are equally
  spaced. Note that the top of an aquifer is also the bottom of a shale, and
  the same location is not entered twice. In other words, with the default
  values for ``NumZPointsWithinAquifers`` (10) and ``NumZPointsWithinShales`` (3)
  a z-grid will have ten points from the bottom to the top of an aquifer, then a
  point in the middle of the overlying shale (point 2 of 3 across the shale),
  and then ten points from the bottom to the top of the overlying aquifer
  (etc.). In this example, including points 1 and 3 for the shale would be
  redundant because those points are included for the aquifers below and above
  the shale.

Below, we show two examples of ``TTFD`` plots setup in the ``Plots``
section of a *.yaml* file. The first plot (*pH_Minimum_Input*) has only the entries
required to set up the ``TTFD`` plot type: ``PlumeType`` and ``ComponentNameList``.
The second plot (*TDS_All_Options_Specified.tiff*) includes all optional entries
for the TTFD plot type. Although there are only two plot entries included,
each entry can result in the creation of multiple figures (e.g., earliest plume timings,
TTFD from monitoring locations, and plume probabilities for each model realization). Note
that all entries for the ``TTFD`` plot type are indented under ``TTFD`` which is
indented under the figure name.

.. code-block:: python
   :lineno-start: 1

    Plots:
        pH_Minimum_Input:
            TTFD:
                PlumeType: pH
                ComponentNameList: [FutureGen2AZMI1, FutureGen2Aquifer1]
        TDS_All_Options_Specified.tiff:
            TTFD:
                PlumeType: TDS
                ComponentNameList: [FutureGen2AZMI1, FutureGen2Aquifer1]
                FigureDPI: 300
                MonitoringLocations:
                    coordx: [100, 200]
                    coordy: [100, 200]
                    coordz: [-407.5, -407.5]
                    HorizontalWindow: 1
                    VerticalWindow: 5
                PlotInjectionSites: True
                InjectionCoordx: 50
                InjectionCoordy: 50
                SpecifyXandYLims:
                    xLims: [-25, 700]
                    yLims: [-25, 700]
                NumZPointsWithinAquifers: 10
                NumZPointsWithinShales: 3
                xGridSpacing: 5
                yGridSpacing: 5
                SpecifyXandYGridLims:
                    gridXLims: [25, 650]
                    gridYLims: [25, 650]
                WriteDreamOutput: False
                SaveCSVFiles: True

For examples of TTFD plots, see *ControlFile_ex39.yaml* to *ControlFile_ex43.yaml*.

GriddedRadialMetric
~~~~~~~~~~~~~~~~~~~
The ``GriddedRadialMetric`` plot type produces map view images of a gridded
radial metric. The ``GenericAquifer`` produces four kinds of gridded
radial metrics: **r_coordinate**, **z_coordinate**, **Dissolved_CO2_mass_fraction**,
and **Dissolved_salt_mass_fraction**. Regions of an aquifer with dissolved |CO2| and
dissolved salt mass fractions exceeding the corresponding mass fraction threshold
are included in the plume volumes for the corresponding plume type. Those plume
volumes can be visualized with the ``TTFD`` plot type. The ``GriddedRadialMetric``
plot type, however, can show more general changes in dissolved |CO2| and salt
mass fractions (e.g., seeing changes in mass fractions below the plume definition
thresholds).

The ``GriddedRadialMetric`` plot type has three required entries: ``ComponentNameList``,
``ZList``, and ``MetricName``. ``ComponentNameList`` is discussed above.

* ``ZList`` - the depths (m) at which to evaluate the radial metric output. The depth in the
  radial grid (e.g., **z_coordinate** from ``GenericAquifer``) that is closest to each value
  entered will be used. String inputs representing the bottom of a unit can also be
  provided. For example, the bottom and top depths of aquifer 2 can be set up by entering
  ``ZList: [aquifer2Depth, shale3Depth]``. Shale 3 is on top of aquifer 2, so the bottom
  depth of shale 3 is the top depth of aquifer 2. Note that numeric values given for
  ``ZList`` (not string inputs like shale2Depth) are taken as being negative when they
  represent a depth beneath the surface (e.g., ``ZList: [-500, -400]`` for depths of 500 m
  and 400 m).

* ``MetricName`` - the name of the gridded radial metric to plot. For example, the
  ``GenericAquifer`` outputs **Dissolved_CO2_mass_fraction** or **Dissolved_salt_mass_fraction**
  can be provided for ``MetricName``. When plotting these metrics, the ``GenericAquifer``
  component used must also produce **r_coordinate** and **z_coordinate** outputs.

The ``GriddedRadialMetric`` plot type has 11 optional entries: ``MinValue``,
``DegreeInterval``, ``Realization``, ``TimeList``, ``PlotInjectionSites``, ``InjectionCoordx``,
``InjectionCoordy``, ``SpecifyXandYLims``, ``SaveCSVFiles``, ``EqualAxes``, and ``FigureDPI``.
All of these entries except for ``MinValue``, ``DegreeInterval``, and ``EqualAxes`` are
discussed above.

* ``MinValue`` - the minimum value used for the colorbar on the figures. Any values beneath
  this minimum will not be displayed graphically, but the entire range of values is still
  displayed in the title of each figure. This parameter has a significant impact on
  ``GriddedRadialMetric`` figures. For example, the **Dissolved_CO2_mass_fraction** and
  **Dissolved_salt_mass_fraction** outputs saved by a ``GenericAquifer`` for a time of 0 years
  will all have values of zero. The outputs saved at other times, however, can have very low but
  nonzero values. The **Dissolved_CO2_mass_fraction** values can be as low as 5.0e-3 at the highest
  radii, while the **Dissolved_salt_mass_fraction** values can be as low as 1.0e-30. If the
  ``MinValue`` provided is zero, then the figures created will be zoomed out to encompass such
  low values at the highest radii evaluated (about 77.5 km). These large extents will make the
  figures visually unclear. For these figures to be useful, one should specify a ``MinValue``
  that is high enough to enable the figure to focus on the area of interest (i.e., near the
  component's location) but low enough to not exclude too much of the output data. We recommend
  using a ``MinValue`` of 0.002 when evaluating **Dissolved_salt_mass_fraction** (10 times lower than
  the default **dissolved_salt_threshold** of 0.02) and 0.01 when evaluating **Dissolved_CO2_mass_fraction**
  (equal to the default **dissolved_salt_threshold** of 0.01). If ``MinValue`` is not entered, these
  values will be used as the defaults for the corresponding output type (dissolved salt or
  dissolved |CO2|). If all of the times evaluated only have values less or equal to ``MinValue``,
  then one figure will be made. This figure has a title that includes 'All Model Times.' Note
  that the .csv files saved when ``SaveCSVFiles`` is set to ``True`` will only include
  values above ``MinValue``.

* ``DegreeInterval`` - the interval (degrees) used to create a map-view image from the radial
  output. The accepted values are 1, 5, 10, 15, 30, and 45. If ``DegreeInterval`` is not entered,
  the default values is 15 degrees.

* ``EqualAxes`` - the option to force the x and y axes to cover the same distances (for an equal
  aspect ratio). The acceptable values are ``True`` or ``False``, and the default value is ``True``.

Note that although the ``Realization`` entry for the ``GriddedRadialMetric`` plot type
follows the indexing conventions of Python (i.e., ``Realization: 0`` for the first realization),
the figure files and .csv files saved by the ``GriddedRadialMetric`` plot type will present the
simulation number as ranging from one to the total number of realizations (e.g., Simulation 1
instead of Simulation 0).

Below, we show two examples of ``GriddedRadialMetric`` plot entries in a control file.
The first entry (*Min_Input_Dissolved_Salt*) uses the minimum input required for the
``GriddedRadialMetric`` plot type. The second entry (*All_Input_Dissolved_Salt*)
uses all entries available for the ``GriddedRadialMetric`` plot type.

.. code-block:: python
   :lineno-start: 1

    Plots:
        Min_Input_Dissolved_Salt:
            GriddedRadialMetric:
                ComponentNameList: [GenericAquifer1]
                MetricName: Dissolved_salt_mass_fraction
                ZList: [aquifer2Depth]
        All_Input_Dissolved_Salt:
            GriddedRadialMetric:
                ComponentNameList: [GenericAquifer1]
                MetricName: Dissolved_salt_mass_fraction
                ZList: [aquifer2Depth, shale3Depth]
                TimeList: [1, 5, 10, 15, 20]
                MinValue: 0.002
                FigureDPI: 300
                PlotInjectionSites: True
                InjectionCoordx: 100
                InjectionCoordy: 100
                DegreeInterval: 1
                Realization: 0
                EqualAxes: True
                SaveCSVFiles: True
                SpecifyXandYLims:
                    xLims: [-200, 400]
                    yLims: [-200, 400]

For examples of ``GriddedRadialMetric`` plots, see *ControlFile_ex53a.yaml* to
*ControlFile_ex53d.yaml*

AtmPlumeSingle
~~~~~~~~~~~~~~
The ``AtmPlumeSingle`` plot type produces map view images depicting how |CO2| leakage
at the surface creates atmospheric |CO2| plumes. These images are created for each
time step during one realization of a simulation. Note that simulations using the
Latin Hypercube Sampling (``lhs``) or Parameter Study (``parstudy``) analysis types have
many realizations, while a simulation using a ``forward`` analysis type only has one
realization. For ``AtmPlumeSingle`` plot type with ``lhs`` or ``parstudy``
simulations, the visualization corresponding to the realization of interest
can be setup with the ``Realization`` entry in the *.yaml* file (discussed above).
Note that using the ``AmtPlumeSingle`` plot type requires the use of an AtmosphericROM
component.

Here is an example of the ModelParams section from *ControlFile_ex40.yaml*, where the
number of LHS realizations is set as ``siz: 30``.

.. code-block:: python
   :lineno-start: 1

    ModelParams:
        EndTime: 15.
        TimeStep: 1
        Analysis:
            type: lhs
            siz: 30
        Components: [LookupTableReservoir1, MultisegmentedWellbore1,
                     FutureGen2AZMI1, FutureGen2AZMI2]
        OutputDirectory: output/output_ex40_{datetime}
        Logging: Info

The produced figures show the source of the |CO2| leak as a red circle and the plume
as a blue circle. The source location(s) are set by the x and y coordinate(s) of
the component that the ``AtmosphericROM`` is connected to. For example, in
*ControlFile_ex9a.yaml*, the ``AtmosphericROM`` component is connected
to an ``OpenWellbore`` component and the ``OpenWellbore`` component has
its locations entered with ``coordx`` and ``coordy``. The ``coordx`` and
``coordy`` values serve as the coordinates of sources for the ``AtmosphericROM``
component. In the ``AtmPlumeSingle`` figures, the ``coordx`` and ``coordy``
values are shown as the |CO2| sources. In the final figures the plumes are labeled
as *Critical Areas* because the area is defined as being within the **critical_distance**
output (from an ``AtmosphericROM``) from the corresponding source. The critical areas are,
therefore, the areas in which the |CO2| concentrations exceed the value defined
by the parameter **CO_critical**. The **critical_distance** is the radius of each plume
circle shown in ``AtmPlumeSingle`` plots, and this **critical_distance** is also
displayed on the figure with text.

Note that when multiple atmospheric plumes overlap enough, they will be displayed
as one plume. The source shown will be between the sources of each individual
plume.

``AtmosphericROM`` components can be provided with receptor locations, which are meant to
represent home or business locations where people will be present. If receptors
are provided and the *.yaml* input for the ``AtmPlumeSingle`` includes the entry
``PlotReceptors: True``, then receptor locations will be shown.

The ``AtmPlumeSingle`` plot type can have the following optional entries ``Realization``,
``PlotReceptors``, ``PlotInjectionSites``, ``InjectionCoordx``, ``InjectionCoordy``,
``SpecifyXandYLims``, and ``FigureDPI``. All of these entries except for ``PlotReceptors``
are described above.

* ``PlotReceptors`` - option to plot receptor locations (default is ``False``). The
  only acceptable values are ``True`` or ``False``. If the receptors are far away from
  the source location(s) and/or the injection site, plotting the receptors may
  cause the x and y limits to be spread too far. The plumes may then be
  difficult to see.

Below is an example of the ``AtmPlumeSingle`` plot input in a *.yaml* control file.
Note that ``InjectionCoordx`` and ``InjectionCoordy`` only have to be provided when
using a ``LookupTableReservoir`` and setting ``PlotInjectionSites: True``.

.. code-block:: python
   :lineno-start: 1

    Plots:
        ATM_single:
            AtmPlumeSingle:
                Realization: 10
                FigureDPI: 300
                PlotInjectionSites: True
                InjectionCoordx: 3.68e4
                InjectionCoordy: 4.83e4
                PlotReceptors: True
                SpecifyXandYLims:
                    xLims: [3.58e4, 3.78e4]
                    yLims: [4.73e4, 4.93e4]

For examples of ``AmtPlumeSingle`` plots, see *ControlFile_ex9a.yaml* to
*ControlFile_ex9c.yaml*.

AtmPlumeEnsemble
~~~~~~~~~~~~~~~~
The ``AtmPlumeEnsemble`` plot type can only be used in simulations with Latin
Hypercube Sampling (``lhs``) or Parameter Study (``parstudy``) analysis types. This
plot type involves concepts similar to those as those of the ``AtmPlumeSingle``
plot type. While the ``AtmPlumeSingle`` plot type dislays the critical areas for
one realization, the ``AtmPlumeEnsemble`` plot type displays the probability of
critical areas occuring in the domain. These probabilities are calculated with
the results from all realizations of the ``lhs`` or ``parstudy`` simulation. The
probabilities specifically represent the likelihood of |CO2| plume concentrations
exceeding the threshold set with the **CO_critical** parameter for ``AtmosphericROM``
components. The probabilities are shown as gridded data. The ``AtmPlumeEnsemble``
plot type is available only when example setup includes an ``AtmosphericROM`` component.

The ``AtmPlumeEnsemble`` plot type has the optional entries ``PlotReceptors``,
``PlotInjectionSites``, ``InjectionCoordx``, ``InjectionCoordy``, ``SpecifyXandYGridLims``,
``xGridSpacing``, ``yGridSpacing``, ``SpecifyXandYLims``, and ``FigureDPI``. All of these
entries are described above.

Below is an example of a ``AtmPlumeEnsemble`` plot entry in a *.yaml* file:

.. code-block:: python
   :lineno-start: 1

    Plots:
        ATM_Ensemble.tiff:
            AtmPlumeEnsemble:
                FigureDPI: 300
                PlotInjectionSites: True
                InjectionCoordx: 200
                InjectionCoordy: 200
                PlotReceptors: False
                xGridSpacing: 1
                yGridSpacing: 1
                SpecifyXandYGridLims:
                    gridXLims: [-100, 300]
                    gridYLims: [-100, 300]
                SpecifyXandYLims:
                    xLims: [-125, 325]
                    yLims: [-125, 325]

For examples of ``AmtPlumeEnsemble`` plots, see *ControlFile_ex9a.yaml* and
*ControlFile_ex9c.yaml*.
