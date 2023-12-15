.. _workflow:

Workflows in NRAP-Open-IAM
==========================

NRAP-Open-IAM simulations in the control file interface or GUI can be set up manually, with the user 
handling all details for the selected components (e.g., parameters, component connections, and outputs). 
Alternatively, the user can set up a simulation by selecting a workflow type. Each workflow type 
is made for a specific type of analysis. An analysis type can require certain types of components, 
and a workflow will automatically add and configure the required components. Additionally, the workflow 
will automatically set up the plots required for the analysis. By building a simulation from the blueprint 
of a workflow, the system model does not have to start from a blank slate. This approach is meant to 
streamline the analyses that are frequently conducted for geologic carbon storage sites.

For example, the system model created by the ``AoR`` workflow will always include a reservoir component, 
a wellbore component, and an aquifer component because those component types are required for an ``AoR`` 
analysis. The user does not need to specify that the wellbore component is connected to the reservoir 
component (i.e., it uses pressures and |CO2| saturations from the reservoir component), because that 
arrangement is assumed within the ``AoR`` workflow. By automatically handling such considerations, the 
``AoR`` workflow allows ``AoR`` analyses to be conducted more quickly and easily.

Workflow types
--------------

The workflows currently available are ``LeakageAssessment``, ``AoR``, and ``TTFD``. We discuss each of 
these workflow types below.

LeakageAssessment Workflow
--------------------------

The ``LeakageAssessment`` workflow is used to assess |CO2| and brine leakage rates given output from a 
reservoir component and the presence of (a) wellbore(s). For example, the wellbore component used could 
represent a legacy well. This workflow automatically adds and sets up a reservoir component and a 
wellbore component - no aquifer component is included. The ``LeakageAssessment`` workflow then makes 
four time series plots. Two of these plots show reservoir pressure and |CO2| saturation at the wellbore 
location(s) given. The other two plots show the leakage rates of |CO2| and brine through the wellbore(s) 
and to specific aqufier. Although this system model configuration is relatively simple, the 
``LeakageAssessment`` workflow is meant to provide convenience while also preventing the errors that can 
occur when manually setting up a system model (e.g., mistyping a component's name).

AoR Workflow
------------

The ``AoR`` workflow conducts an area of review analysis for a geologic carbon storage site. An ``AoR`` 
analysis evaluates the maximum potential impacts on the reservoir and an aquifer. The reservoir impacts 
considered are pressures and |CO2| saturations. The aquifer impacts considered are the contaminant plume 
volumes (e.g., pH and TDS plumes) that could arise if a leakage pathway was located at a particular 
location. Many hypothetical leakage pathways (wells) are considered, with their positions distributed 
across the entire study area. The ``AoR`` is determined as the area in which the hypothetical wells 
could cause the aquifer to be contaminated. In other words, the ``AoR`` includes any well location that, 
in the simulation, had a nonzero plume volume in the aquifer, a nonzero |CO2| saturation in the 
reservoir, or a reservoir pressure exceeding the critical pressure. The critical pressure is 
calculated as the pressure required to drive flow from the reservoir into the aquifer being 
considered. A critical pressure can be automatically calculated based on the reservoir depth, 
aquifer depth, and a specified brine density, or it can be given as a specific value.

Note that the ``AoR`` workflow is not the same as the ``AoR`` plot type (discussed in section
:ref:`cfi_visualization`). While the ``AoR`` plot type only shows the results for one metric (e.g.,
reservoir pressure), the ``AoR`` workflow assesses all four metrics considered (four separate ``AoR``
plots) and creates an overall ``AoR`` that reflects all of these metrics.

``AoR`` analyses can be conducted with any type of wellbore component, but we recommend using the 
``OpenWellbore`` component. This component portrays an uncemented wellbore, which represents a 
worst-case scenario for leakage risks. When an ``AoR`` might contain undocumented legacy wells, it can 
be best to assume a conservatively high portrayal of leakage risks. For example, if the undocumented 
legacy wells were assumed to be cemented, but some are uncemented or some have a higher effective 
permeability than assumed in the simulations, then the leakage risks would be underestimated and 
the resulting ``AoR`` would be too small. Additionally, the ``OpenWellbore`` component can utilize a 
critical pressure in leakage calculations. Other wellbore components (e.g., ``MultisegmentedWellbore`` 
and ``CementedWellbore``) do not use a critical pressure. When using one of these wellbore components, 
the critical pressure would still be shown in the pressure ``AoR`` plot (i.e., one of the five map-view 
.figures created by the ``AoR`` workflow - this figure only shows results from the reservoir component). 
The critical pressure would not, however, be reflected in the pH and TDS plume volume ``AoR`` plots (two 
of the five map-view figures). In other words, these plume volumes result from |CO2| and brine 
leakage into the aquifer, and this leakage is calculated by an ``OpenWellbore``, 
``MultisegmentedWellbore``, or ``CementedWellbore`` component. The ``MultisegmentedWellbore`` and 
``CementedWellbore`` components will calculate leakage in their default manner, without a 
user-specified critical pressure. Nonetheless, the user may want to use one of these components 
if they have constraints on the effective permeabilities of legacy wells in the area.

If the critical pressure (|Pcrit|) is not given as a specific value, it is calculated with the approach 
shown in section :ref:`equations`. When using the ``AoR`` workflow, the critical pressure input provided 
is used in both (1) the map-view figure of reservoir pressures, which impacts the overall ``AoR`` figure, 
and (2) the leakage calculations of ``OpenWellbore`` components (if one is used). If the ``AoR`` workflow
is using a wellbore component with a brine density parameter (e.g., **brineDensity** parameter of the
``OpenWellbore`` and ``MultisegmentedWellbore`` components), then that parameter will be used for the
brine density value in the critical pressure equation. If the wellbore component does not have a brine density
parameter (e.g., ``CementedWellbore``), then a brine density can be specified with the ``BrineDensity``
entry (discussed further below).

While the ``AoR`` plot type can be used with a ``TheisReservoir``, the ``AoR`` workflow cannot be used
with a ``TheisReservoir``. When using a ``TheisReservoir``, ``AoR`` plots can only be made for the 
**pressure** output (see *ControlFile_ex47*). The ``AoR`` workflow is meant to incorporate data from all
of the metrics that can be used with the ``AoR`` plot type (i.e., plume volumes and **CO2saturation**),
so the ``AoR`` workflow cannot be used with a ``TheisReservoir``.

TTFD Workflow
-------------

The ``TTFD`` workflow conducts a time to first detection analysis. This analysis type examines 
the spread of a contaminant plume through an aquifer, with multiple map-view figures showing 
how the plume spreads over time. The time to first detection is obtained by providing the 
locations and depths of monitoring wells. The time at which a plume first reaches a monitoring 
well is taken as the ``TTFD``. If the monitoring well locations used do not detect a plume, then 
the locations should be changed. This workflow can also produce the files used as input for 
DREAM (Designs for Risk Evaluation and Management), a tool designed to optimize monitoring network 
design at geological carbon storage sites.

The ``TTFD`` plot type, which creates map-view images of aquifer plume timings across the domain, 
exists separately from the ``TTFD`` workflow. The ``TTFD`` workflow is meant to offer convenience 
by automatically setting up the components, component connections, and plot entries (in a 
*.yaml* file) required for use of the ``TTFD`` plot type.

Setup of a Workflow in the Control File Interface
-------------------------------------------------

In a standard control file, there are separate sections for stratigraphy, model parameters, 
each component, and the plots to be produced. A control file using a workflow, however, can 
accomplish complex analyses while only using three sections: a stratigraphy section, a model 
parameters section, and a workflow section. We show an example of this setup below.

.. code-block:: python
   :lineno-start: 1

    #-------------------------------------------------
    ModelParams:
        EndTime: 20
        TimeStep: 1.0
        Analysis: forward
        Components: []
        OutputDirectory: 'output/output_ex55a_{datetime}'
        GenerateOutputFiles: False
        GenerateCombOutputFile: False
        Logging: Debug
    #-------------------------------------------------
    Stratigraphy:
        numberOfShaleLayers:
            value: 5
            vary: False
        shale1Thickness:
            value: 198.7
            vary: False
        shale2Thickness:
            value: 74.4
            vary: False
        shale3Thickness:
            value: 110.3
            vary: False
        aquifer1Thickness:
            value: 33.2
            vary: False
        aquifer2Thickness:
            value: 84.1
            vary: False
        reservoirThickness:
            value: 7.0
            vary: False
    #-------------------------------------------------
    Workflow:
        Type: LeakageAssessment
        Options:
            ReservoirComponentType: AnalyticalReservoir
            WellboreComponentType: OpenWellbore
            ReservoirOptions:
                InjectionWell:
                    coordx: 10
                    coordy: 50
                Parameters:
                    injRate: 0.1
            WellboreOptions:
                Locations:
                    coordx: [500]
                    coordy: [500]
                Parameters:
                    wellRadius: 0.05
            AquiferName: aquifer2
            FigureDPI: 300
    #-------------------------------------------------

The ``ModelParams`` and ``Stratigraphy`` sections are set up in the same manner as a normal control
file (see section :ref:`control_file`). The ``Workflow`` section contains two entries: ``Type`` and ``Options``. The 
``Type`` entry is followed by one of the available workflow types. Here, the ``Type`` is ``LeakageAssessment``. 
The ``Options`` entry is a dictionary containing more entries within it (i.e., in a *.yaml* control 
file, the contained entries are on a line beneath ``Options`` and preceeded by four additional spaces). 
The remaining details regarding the simulation are entered under ``Options`` (e.g., well locations, 
component parameters, and plotting options). Any options that are not given will be set at their default 
settings. For the sake of brevity, the example above does not show all of the entries that can be provided 
under ``Options``. Some of the ``Options`` entries available depend on the workflow type, and the entries 
specific to each workflow are discussed further below.

The ``Options`` entries that are applicable across multiple workflow types are:

* ``ReservoirComponentType`` - the type of reservoir component to use. Acceptable inputs include 
  ``LookupTableReservoir`` and ``AnalyticalReservoir``. Any analysis conducted for a real site 
  should use a ``LookupTableReservoir`` to incorporate results from high-fidelity reservoir 
  simulations; the remaining reservoir component types are largely intended for testing purposes.

* ``WellboreComponentType`` - the type of wellbore component to use. Acceptable inputs include 
  ``OpenWellbore``, ``MultisegmentedWellbore``, and ``CememntedWellbore``.

* ``AquiferComponentType`` - the type of aquifer component to use. Acceptable inputs include 
  ``GenericAquifer``, ``FutureGen2Aquifer``, ``FutureGen2AZMI``, and ``DeepAlluviumAquifer``.

* ``ReservoirOptions`` - an entry containing other entries related to the reservoir component. All 
  of the input provided in the reservoir component's section of a standard control file are included 
  here (e.g., the ``Parameters`` and ``InjectionWell`` entries).

* ``WellboreOptions`` - an entry containing other entries related to the wellbore component. All 
  of the input provided in the wellbore component's section of a standard control file are included 
  here (e.g., the ``Locations``, ``LeakTo``, ``Controls``, ``ThiefZone``, and ``Parameters`` entries).
  For examples of how to set up entries like ``Locations`` and ``Parameters``, see section :ref:`control_file`.
  When using an ``OpenWellbore`` component in a workflow, the default approach is to include the ``Controls`` 
  section with ``critPressureApproach`` set to ``True`` and ``enforceCritPressure`` set to ``False``. 
  With these settings, the ``OpenWellbore`` component will automatically use a critical pressure 
  calculated with the **brineDensity** parameter, the bottom depth of the aquifer receiving leakage, 
  and the top depth of the reservoir.

* ``AquiferOptions`` - an entry containing other entries related to the aquifer component. All 
  of the input provided in the aquifer component's section of a standard control file are included 
  here (e.g., the ``Parameters`` entry).

* ``AquiferName`` - the name of the aquifer being considered in the analysis, formatted as ``aquifer2`` 
  (where ``2`` is replaced by the desired aquifer number). If this entry is provided under the 
  ``AquiferOptions`` entry, then that input will overwrite any input given for ``AquiferName`` under 
  the ``Options`` entry. If not provided, the default setting is the highest aquifer.

* ``PlotInjectionSites`` - the option to plot injection sites on any map-view figures created. The only 
  acceptable values are ``True`` or ``False``. The default setting is ``False``.

* ``InjectionCoordx`` - value or list of values for the x coordinate(s), or easting distance(s), of the 
  injection site(s) (default is None). The values are given in meters. This entry must be provided when 
  using a ``LookupTableReservoir``, as that component type does not have an *.injX* attribute. Other 
  reservoir types like ``AnalyticalReservoir`` can be displayed without an ``InjectionCoordx`` entry.

* ``InjectionCoordy`` - value or list of values for the y coordinate(s), or northing distance(s), of the 
  injection site(s) (default is None). The values are given in meters. This entry must be provided when 
  using a ``LookupTableReservoir``, as that component type does not have an *.injY* attribute. Other 
  reservoir types like ``AnalyticalReservoir`` can be displayed without an ``InjectionCoordy`` entry.

* ``FigureDPI`` - the dots-per-inch (DPI) of the figure(s) produced (default is 100). Larger DPIs 
  will create high-resolution, high-quality figures, but the file sizes are also larger. File size 
  and quality are also influenced by the extension used (e.g., *.png* or *.tiff*). Recommended 
  ``FigureDPI`` values are between 100 and 300.

* ``FigureSize`` - the width and height of the figures produced in inches as a list of length two 
  (e.g., ``FigureSize: [14, 8]```), where the first number is the width and the second number is the height.

* ``FigureName`` - the name of the figure to be produced. If the name is given with an extension, 
  then the figure will be saved using that file extension (e.g., ``FigureName: Overall_AoR_Plot.tiff`` 
  for a .tiff file named Overall_AoR_Plot). If this entry is not provided, the default figure names 
  will be used. If this entry is provided for a workflow that produces multiple figures (e.g., the 
  ``LeakageAssessment`` and ``TTFD`` workflows), the name will not be used but any extension included 
  will be used.

* ``AutomateResCompSetup`` - the option specifying whether the workflow should automatically set up 
  the section for the reservoir component in the *.yaml* file. The only acceptable values are ``True`` 
  or ``False``, and the default setting is ``True``.

* ``AutomateWellCompSetup`` - the option specifying whether the workflow should automatically set up 
  the section for the wellbore component in the *.yaml* file. The only acceptable values are ``True`` 
  or ``False``, and the default setting is ``True``.

* ``AutomateAqCompSetup`` - the option specifying whether the workflow should automatically set up 
  the section for the aquifer component in the *.yaml* file. The only acceptable values are ``True`` 
  or ``False``, and the default setting is ``True``.

* ``AutomatePlotsSetup`` - the option specifying whether the workflow should automatically set up 
  the Plots section of the *.yaml* file. The only acceptable values are ``True`` or ``False``, and 
  the default setting is ``True``. If set to ``False``, then the ``Plots`` section of the *.yaml* 
  control file must be set up by the user.

Setup of the LeakageAssessment Workflow in the Control File Interface
---------------------------------------------------------------------

To use the ``LeakageAssessment`` workflow, set the ``Type`` entry of the ``Workflow`` section to 
``LeakageAssessment``.

The ``LeakageAssessment`` workflow requires the ``Locations`` input under ``WellboreOptions``:

.. code-block:: python
   :lineno-start: 1
    
    #-------------------------------------------------
    Workflow:
        Type: LeakageAssessment
        Options:
            PlotType: TimeSeriesStats
            FigureDPI: 300
            UseMarkers: True
            Subplot:
                Use: True
                NumCols: 2
            FigureSize: [14, 8]
            FigureName: Name.tiff   # Only the .tiff extension from this will be used
            ReservoirComponentType: AnalyticalReservoir
            ReservoirOptions:
                InjectionWell:
                    coordx: 0
                    coordy: 0
                Parameters:
                    injRate: 0.5
                    reservoirRadius: 2000
                    logResPerm: -13.5
                    brineDensity: 1100
            WellboreComponentType: OpenWellbore
            WellboreOptions:
                Locations:
                    coordx: [100, 200]
                    coordy: [100, 200]

Here, the ``Locations`` entry is given ``coordx`` and coordy`` lists. Any valid input for ``Locations`` 
can be given under ``WellboreOptions`` (e.g., the ``file``, ``grid``, or ``range`` options; see 
*ControlFile_ex30.yaml*).

In addition to the optional entries that are applicable to multiple workflows (discussed above), the 
``LeakageAssessment`` workflow has the optional entries ``PlotType``, ``Subplot``, ``UseMarkers``, 
``UseLines``, ``VaryLineStyles``, and ``FigureSize`` (all of which are entered under ``Options``). 
All of these entries except for ``PlotType`` are options for the ``TimeSeries`` plot type that control 
the appearance of a time series plot (e.g., specifying whether to use multiple subplots). For more 
information about those entries, see section :ref:`cfi_visualization`.

* ``PlotType`` - option specifying the type of time series plot. The options are ``TimeSeries``, 
  ``TimeSeriesStats``, and ``TimeSeriesAndStats``. For more details regarding these plot types, see 
  section :ref:`cfi_visualization`.

The figures produced will focus on the leakage received by the aquifer specified through the 
``AquiferName`` entry discussed above. If that entry is not provided, then the default setting is the 
highest aquifer.

For examples of the ``LeakageAssessment`` workflow in the control file interface, see *ControlFile_ex58a.yaml* 
to *ControlFile_ex58c.yaml*.

Setup the AoR Workflow in the Control File Interface
----------------------------------------------------

To use the ``AoR`` workflow, set the ``Type`` entry of the ``Workflow `` section to ``AoR``.

When using the ``AoR`` workflow, we recommend setting ``GenerateOutputFiles`` and ``GenerateCombOutputFile`` 
to ``False`` in the ``ModelParams`` section of the *.yaml* file. The large number of wellbore locations commonly 
used for ``AoR`` plots causes a large number of output files. A reservoir and aquifer component is created for
each wellbore location, and every component will have its output saved. The ``.csv`` files created by the ``AoR`` 
workflow contain all of the necessary information and these files are much smaller in size.

In addition to the optional entries that are applicable to multiple workflows (discussed above), the 
``AoR`` Workflow has the optional entries  ``TimeList``, ``CriticalPressureMPa``, and ``BrineDensity``. These 
entries can be provided under ``Options``. The ``TimeList`` and ``BrineDensity`` entries have the same impact on 
both the ``AoR`` plot type and the ``AoR`` workflow; for more details on these entries, refer to the section
:ref:`cfi_visualization`. The ``CriticalPressureMPa`` entry has a different impact when used for an ``AoR``
plot entry or for the ``AoR`` workflow, however, so we describe the ``CriticalPressureMPa`` entry for the
``AoR`` workflow below.

* ``CriticalPressureMPa`` - this entry controls how critical pressure is handled in the ``AoR`` analysis. 
  The entry can either be given as ``Calculated`` or as a number representing a critical pressure in 
  MPa (e.g., ``CriticalPressureMPa: 20.5`` for 20.5 MPa). The ``AoR`` plot type also has an entry called 
  ``CriticalPressureMPa``, but when that entry is used for an ``AoR`` plot it only has an impact if the ``AoR`` 
  plot analyzes the *pressure* metric (i.e., no impact on plots evaluating |CO2| saturations or aquifer plume 
  volumes). When ``CriticalPressureMPa`` is given for the ``AoR`` workflow, however, this setting impacts both 
  the critical pressures used by an ``OpenWellbore`` component (if one is used) and the critical pressures used 
  in ``AoR`` plots evaluating the *pressure* metric.

Below, we show an example of an ``AoR`` workflow in the control file interface. Only the ``Workflow`` 
section is shown (the ``ModelParams`` and ``Stratigraphy`` sections are not shown).

.. code-block:: python
   :lineno-start: 1
    
    #-------------------------------------------------
    Workflow:
        Type: AoR
        Options:
        FigureName: Overall_AoR_Plot.tiff
        CriticalPressureMPa: Calculated
        ReservoirComponentType: LookupTableReservoir
        ReservoirOptions:
            FileDirectory: source/components/reservoir/lookuptables/FutureGen2/1008_sims
            TimeFile: time_points.csv
            ParameterFilename: parameters_and_filenames_trunc.csv
            Interpolation2D: False
            Parameters:
                index: 1
        WellboreComponentType: OpenWellbore
        WellboreOptions:
            Locations:
                file: source/components/reservoir/lookuptables/FutureGen2/1008_sims/fg1.csv
                read_z_values: True
                thin_point_density: True
                min_x_spacing: 10000
                min_y_spacing: 10000
            Parameters:
                wellRadius: 0.05
                logReservoirTransmissivity: -10.0
                logAquiferTransmissivity: -10.0
                brineSalinity: 0.0475
                brineDensity: 1075
        AquiferName: aquifer4
        AquiferComponentType: FutureGen2Aquifer
        AquiferOptions:
            Parameters:
                por: 0.18
                log_permh: -11.92
                log_aniso: 0.3

Unlike the ``LeakageAssessment`` and ``TTFD`` workflows, the ``AoR`` workflow can be run without 
explicitly setting the locations of wellbores. Specifically, the user does not have to enter x and y 
values through the ``Locations`` entry under ``WellboreOptions``, as discussed above. The user can enter 
that information, but if that input is not provided then there are default settings for the ``AoR`` workflow.
When using a ``LookupTableReservoir``, the default x and y values for the wellbore components are taken 
directly from the .csv files tied to the ``LookupTableReservoir`` component (i.e., by using the ``file`` entry 
under ``Locations``; see *ControlFile_ex30.yaml*). Because the point densities of reservoir models generally 
become quite high near an injection site, the default approach is to thin the point density by enforcing a 
minimum point spacing of 20 km in the x and y directions (the ``min_x_spacing`` and ``min_y_spacing`` entries 
shown above). This large spacing is the default setting because these simulations can have long run times. 
Reducing the point density of wellbore locations reduces the number of reservoir, wellbore, and aquifer 
components used. If the reservoir component is not a ``LookupTableReservoir``, then the default wellbore 
locations are arranged in a six by six grid (total of 36 points) with x and y values ranging from -50 km to 
50 km.

To manually define wellbore locations with a grid, use the ``grid`` entry under ``Locations``. The ``grid`` 
entry contains the entries ``xmin``, ``xmax``, ``xsize``, ``ymin``, ``ymax``, and ``ysize``. The minimum 
x and y values (in meters) are defined by ``xmin`` and ``ymin``, respectively. The maximum x and y values 
(in meters)are defined by ``xmax`` and ``ymax``, respectively. Finally, the number of grid points in the 
x and y directions are defined by ``xsize`` and ``ysize``, respectively.

.. code-block:: python
   :lineno-start: 1
    
        WellboreOptions:
            Locations:
                grid:
                    xmin: -7500
                    xmax: 7500
                    xsize: 14
                    ymin: -7500
                    ymax: 7500
                    ysize: 14

In a normal control file (i.e., one that does not use a workflow), an ``OpenWellbore`` component can 
only use a specific critical pressure (as opposed to a calculated one) if the user provides the **critPressure** 
parameter and sets both the ``critPressureApproach`` and ``enforceCritPressure`` entries under ``Controls`` to 
``True``. When using the ``AoR`` workflow that includes an ``OpenWellbore`` component, setting the 
``CriticalPressureMPa`` entry to a specific value (e.g., ``CriticalPressureMPa: 22.07`` for 22.07 MPa) also 
sets the critical pressure for the ``OpenWelbore`` component. In this case, you do not need to manually set 
the **critPressure** parameter or the entries under ``Controls`` for the ``OpenWellbore`` component. The workflow 
automatically changes these settings, with the ``CriticalPressureMPa`` value being converted to |Pa| for the 
**critPressure** parameter. See *ControlFile_ex56c.yaml* for an example where the critical pressure is set 
to a specific value.

Note that the ``CriticalPressureMPa`` input will not impact the leakage calculations of wellbore components 
like ``MultisegmentedWellbore`` or ``CementedWellbore``, which do not use a critical pressure.

For examples of the ``AoR`` workflow in the control file interface, see *ControlFile_ex55a.yaml* to 
*ControlFile_ex56g.yaml*.

Setup the TTFD Workflow in the Control File Interface
-----------------------------------------------------

To use the ``TTFD`` workflow, set the ``Type`` entry of the ``Workflow ``section to ``TTFD``.

The ``TTFD`` Workflow requires the input of ``Locations`` data under ``WellboreOptions``. Any valid 
input for ``Locations`` can be given under ``WellboreOptions`` (e.g., the ``file``, ``grid``, or ``range`` options; 
see *ControlFile_ex30.yaml*).

The ``PlumeType`` entry plays an important role in the TTFD workflow.

* ``PlumeType`` - the type of plume metric being considered. Acceptable values are *Pressure*, *pH*, *TDS*, 
  *Dissolved_CO2*, *Dissolved_salt*, and *CarbonateAquifer*. The plume type used must be compatible with the 
  type of aquifer component used. For example, the default aquifer component is ``GenericAquifer``, and that 
  component only produces *Dissolved_CO2* and *Dissolved_salt* plumes. Using a plume type of *pH* would cause 
  an error when using a ``GenericAquifer``. The default plume type is *pH* when using a ``FutureGen2Aquifer``, 
  ``FutureGen2AZMI``, ``DeepAlluviumAquifer``, or ``DeepAlluviumAquiferML`` component. The default plume type 
  is *Dissolved_CO2* when using a ``GenericAquifer`` component. And the default plume type is ``CarbonateAquifer`` 
  when using a ``CarbonateAquifer`` component.

Note that certain components define plume extents with fixed thresholds. For example, ``FutureGen2Aquifer`` and 
``FutureGen2AZMI`` components define *pH* plumes as areas with changes in pH of greater than 0.2. In contrast, 
certain components allow the user to set the thresholds for plume definition. For example, the ``GenericAquifer`` 
defines *Dissolved_CO2* plumes as areas with |CO2| mass fractions exceeding the **dissolved_co2_threshold** parameter. 
The ``GenericAquifer`` also allows the user to set the initial salinity of the aquifer, which impacts how easily 
the **dissolved_salt_threshold** parameter is exceeded and how large the *Dissolved_salt* plumes become.

The ``TTFD`` workflow also has all of the optional entries that apply to the ``TTFD`` plot type (``MonitoringLocations``,  
``WriteDreamOutput``, ``SaveCSVFiles``, ``SpecifyXandYLims, ``SpecifyXandYGridLims``, ``xGridSpacing``, ``yGridSpacing``, 
``NumZPointsWithinAquifers``, and ``NumZPointsWithinShales``). For the ``TTFD`` workflow, all of these entries are 
provided under ``Options`` in the ``Workflow`` section of a .yaml control file. For descriptions of all of these options, 
see the description of the ``TTFD`` plot type in section :ref:`cfi_visualization`.

The input provided under ``MonitoringLocations`` controls the positions and depths of the monitoring wells used. If 
``MonitoringLocations`` is not provided, then the ``TTFD`` workflow will evaluate the evolution of contaminant plumes 
but not show any detection times. If ``WriteDreamOutput`` is set to ``True``, then the ``TTFD`` workflow will save 
*.iam* files containing the plume timing results. These *.iam* files can be used as input to the DREAM tool (Design 
for Risk Evaluation and Management), which was also developed by NRAP.

Below, we show an example of the ``Workflow`` section of a *.yaml* control file using the ``TTFD`` workflow.

.. code-block:: python
   :lineno-start: 1
    
    #-------------------------------------------------
    Workflow:
        Type: TTFD
        Options:
            MonitoringLocations:
                coordx: [200, 800, 200, 800]
                coordy: [800, 200, 800, 200]
                coordz: [-1015, -715, -1015, -715]
            WriteDreamOutput: False
            SaveCSVFiles: False
            FigureDPI: 300
            FigureName: Name.tiff               # specifies .tiff extension
            PlumeType: Dissolved_salt
            PlotInjectionSites: True
            ReservoirComponentType: AnalyticalReservoir
            ReservoirOptions:
                InjectionWell:
                    coordx: 500
                    coordy: 500
                Parameters:
                    injRate: 0.1
                    logResPerm: -13
                    reservoirRadius: 5000
            WellboreComponentType: OpenWellbore
            WellboreOptions:
                Controls:
                    critPressureApproach: True  # Both set to True, so the critPressure
                    enforceCritPressure: True   # parameter will be used
                Locations:
                    coordx: [400, 600]
                    coordy: [600, 400]
                Parameters:
                    wellRadius: 0.05
                    critPressure: 2.207e+7      # equivalent to 22.07 MPa
            AquiferName: aquifer2
            AquiferComponentType: GenericAquifer
            AquiferOptions:
                Parameters:
                    por: 0.118
                    log_permh: -13.39
                    log_aniso: 0.3

For examples of the ``TTFD`` workflow in the control file interface, see *ControlFile_ex57a.yaml* to 
*ControlFile_ex57c.yaml*.
