.. highlight:: python
   :linenothreshold: 3

.. _control_file:

Control File Interface
======================

Control files are a method of providing user input to NRAP-Open-IAM. These files
use a YAML format (extension *.yaml*), and they can be opened with a text editor
(e.g., Notepad). Any line in the control file starting with a pound sign (#) is
a comment that will be ignored by the program. The basic format of the control file is the
name of an entry (e.g., parameter name) followed by a colon, a space, and a value.
Some entries are not directly followed by a value. Instead, the colon is the last
character on the line and then other entries are indented beneath the entry.

The use of NRAP-Open-IAM with a control file is called the control file interface (CFI).
In this section, the use of the CFI is demonstrated with excerpts from *ControlFile_ex1a.yaml*:

.. code-block:: python
   :lineno-start: 1

    #-------------------------------------------------
    # NRAP-Open-IAM Control File Example 1a
    #-------------------------------------------------
    ModelParams:
        EndTime: 50
        TimeStep: 1.0
        Analysis: forward
        Components: [AnalyticalReservoir1,
                     CementedWellbore1]
        OutputDirectory: ../../output/output_ex1a_{datetime}
        Logging: Debug
        GenerateOutputFiles: True
        GenerateCombOutputFile: False
        GenerateStatFiles: False

Here, the first three lines are comments that are not read by the code.
The fourth line defines the keyword ``ModelParams`` which describes parameters of
the system model. The subsequent lines contain parameters of ``ModelParams``.
A ``ModelParams`` section is required in all NRAP-Open-IAM control files. The ``EndTime``
keyword defines the ending time for the simulation in years (50 years in
this example). The ``TimeStep`` parameter defines the length of a time step
(1 year in this example). The type of analysis being run is a ``forward``
(deterministic) simulation. Other possible options for ``Analysis`` parameter
are ``lhs`` for Latin Hypercube Sampling analysis and ``parstudy`` for a parameter
study. See section :ref:`conceptual_model_overview` for more details regarding
these analysis types.

The ``Components`` parameter is a required entry that contains a list of component
model names that are defined later in the file. The component list will always begin
with a square bracket '[' followed by each of the component names that make up the system
separated by a comma ',' and closed by a square bracket ']'. The names of the components
are listed in the order they are supposed to be run. This order can be important,
for example, when a wellbore component needs to run after the reservoir component it
is connected to.

The next keyword ``OutputDirectory`` defines a directory for the output to be
written into. The output directory can be appended with a keyword ``{datetime}``.
When a simulation is run, the ``{datetime}`` keyword will be replaced with the date
and time of the simulation run. Note that the ``{datetime}`` keyword is optional: if it is
omitted subsequent runs of the simulation will overwrite past results.
That is, if there is a need to keep all results from re-running an NRAP-Open-IAM case,
the ``{datetime}`` keyword will easily facilitate this; if re-running an
NRAP-Open-IAM case should overwrite previous results, the ``{datetime}`` keyword should be
omitted. In the output folder, NRAP-Open-IAM places a copy of the input file, a folder
("csv_files") containing all outputs from the component models written to ``.csv`` files,
and files for all of the plots created.

The keyword ``Logging`` defines what level of logging information is written out
to the logging files. Options for ``Logging`` levels are ``Debug``, ``Info``, ``Warning``,
``Error``, and ``Critical``. ``Info`` is the default level (if no ``Logging`` keyword is given, the logging
will be set to ``Info``), and it provides valuable information about when parameters go outside
of permitted ranges and when there are problems in the program. ``Debug`` is a good option
if there are problems with the simulation and more information is required to explore the
causes. A logging level of ``Warning`` will limit the information to only messages deemed
more important (warning messages). Choosing logging levels of ``Error`` or ``Critical``
will limit the messages to only more serious problems, with ``Critical`` messages having
the highest importance.

The ``GenerateOutputFiles``, ``GenerateCombOutputFile``, and ``GenerateStatFiles`` control
the types of output files saved. Each of these entries can be set to ``True`` or ``False``;
if any of these entries are not provided, the default setting is ``True``. If ``GenerateOutputFiles``
is ``True``, the simulation will save results to ``.csv`` files. Each component output will be
saved to its own ``.csv`` file. If ``GenerateCombOutputFile`` is ``True``, all simulation results
will be saved to one file. If the analysis type is ``forward``, the file is called
``simulation_results.csv``. If the analysis type is ``lhs`` or ``parstudy``, the file is called
``output_realizations.csv`` and the parameter values used are saved to another file called
``parameter_realizations.csv``. If ``GenerateStatFiles`` is ``True`` and the analysis type is
``lhs`` or ``parstudy``, the simulation will save two files detailing the statistics of the output
(``output_stats.csv``) and the parameter values used (``parameter_stats.csv``). The statistics
shown are the minimum, maximum, and mean values as well as the standard deviation, variance, and
different percentiles (2.5%, 5%, 25%, 50%, 75%, 95%, and 97.5%; in the files, the 2.5% and 97.5% values
are labelled as the 25th per mille and 975th per mille, respectively, as percentiles do not include
decimals). If the analysis type is ``forward``, the statistics files will not be saved.

.. code-block:: python
   :lineno-start: 1

    #-------------------------------------------------
    Stratigraphy:
        numberOfShaleLayers:
            vary: False
            value: 3
        # Thickness is in meters
        shale1Thickness:
            min: 500.0
            max: 550.0
            value: 525.0
        shale2Thickness:
            min: 450.0
            max: 500.0
            value: 475.0
        shale3Thickness:
            vary: False
            value: 11.2
        aquifer1Thickness:
            vary: False
            value: 22.4
        aquifer2Thickness:
            vary: False
            value: 19.2
        reservoirThickness:
            vary: False
            value: 51.2

The next section of the file is the ``Stratigraphy`` section. This section defines any
model parameters related to the stratigraphy of the |CO2| storage site. See the
:ref:`stratigraphy_component` section of this document for a list of all available
parameters for the ``Stratigraphy`` component. Note that while this example uses a 
``Stratigraphy`` component, there are other types of stratigraphy components available
(e.g., ``DippingStratigraphy`` and ``LookupTableStratigraphy``). A control file must have a
section for one of the stratigraphy component types.

Any parameters for the ``Stratigraphy`` component are defined here with either a deterministic
value or a range to vary over. A fixed value of any given parameter can be specified with the
``vary: False`` and ``value: ###`` specification shown here or simply ``parameterName: ###``
(where ``###`` is the value). If the parameter is meant to vary across different realizations,
the ``min`` and ``max`` entries should be provided under the parameter to represent the
minimum and maximum parameter limits, respectively.

The next sections of the input file define every component model contained in the component
model list specified earlier in the control file (under ``ModelParams``). The first component
in the ``Components`` list is ``AnalyticalReservoir1``, and the component settings are defined
in the following section:

.. code-block:: python
   :lineno-start: 37

    #-------------------------------------------------
    # AnalyticalReservoir1 is a user defined name for component;
    # the type AnalyticalReservoir is the ROM model name
    #-------------------------------------------------
    AnalyticalReservoir1:
        Type: AnalyticalReservoir
        InjectionWell:
            coordx: 10
            coordy: 20
        Parameters:
            injRate: 0.1
        Outputs: [pressure,
                  CO2saturation]

The name *AnalyticalReservoir1* can be replaced with any other name defined by user
(e.g., *ARes1*), but the component will only be included in the system model if its
name is in the ``Components`` list described in the previous section ``ModelParams`` 
(the names used must match). The ``Type`` is a keyword that defines the component model to
be used and must match up with one of the component models currently available in NRAP-Open-IAM.
The ``InjectionWell`` entry specifies the location of the injection well. The ``x`` and ``y`` coordinates
of the injection well are set with the ``coordx`` and ``coordy`` values under ``InjectionWell``,
respectively. If these inputs are not provided, the default injection well location for an
``AnalyticalReservoir`` is ``x`` = 0 |m|, ``y`` = 0 |m|. The ``Parameters`` section defines parameters
of the component model. Descriptions of the parameters available for the user to specify can
be found in the :ref:`components_description` chapter of the current documentation. The component
model parameters are specified in the same fashion as the ``Stratigraphy`` parameters shown above.
The ``Outputs`` entry specifies the observations of the component model that will be output
from the simulation. Please refer to the :ref:`components_description` chapter of this document
to see which parameters and outputs are available for user specification in the CFI.

Generally, dynamic (time-varying) input to component models comes from the
output of other connected component models (e.g., the pressure and saturation
as an input to a wellbore leakage model comes from a reservoir model). In some instances,
there may be a need to study a component model without another attached component
models feeding the input. In this case dynamic input can be specified with
the ``DynamicParameters`` keyword. Under the ``DynamicParameters`` section each
input name is specified followed by a list of values (enclosed in square brackets, "[]")
of the same length as the number of time points (a value for each time point, including
an initial value). See files *ControlFile_ex7a.yaml* and *ControlFile_ex7b.yaml*
for example of control files utilizing dynamic input for some components.

The next section of the input file is similar to the previous section and defines
the next component model *CementedWellbore1*.

.. code-block:: python
   :lineno-start: 47

    #-------------------------------------------------
    CementedWellbore1:
        Type: CementedWellbore
        Connection: AnalyticalReservoir11
        Number: 4
        Locations:
            coordx: [100, 540]
            coordy: [100, 630]
        RandomLocDomain:
            xmin: 150
            xmax: 250
            ymin: 200
            ymax: 300
        Parameters:
            logWellPerm:
                min: -14.0
                max: -12.0
                value: -13.0
        Outputs: [CO2_aquifer1,
                  CO2_aquifer2,
                  CO2_atm,
                  brine_aquifer1,
                  brine_aquifer2]

This part of the example sets up the ``CementedWellbore`` component model named
``CementedWellbore1``. There are four wellbores of this type being added with
``Number: 4``: two of the locations are given in the ``Locations`` part and the
other two are generated randomly within the domain specified in the ``RandomLocDomain``
entry. All coordinate systems are assumed to have units of meters.

Known wellbore coordinates are entered as a comma separated list (i.e., contained in
brackets, "[]"). The x and y coordinates are entered with the ``coordx`` and ``coordy``
entries, respectively. There must be a comma between each coordinate (e.g.,
``coordx: [100, 540]``).

Unknown wellbore locations can be generated by specifying more wellbores (with
``Number``) than the number of known wellbore locations. To control the locations
of randomly placed wells, a ``RandomLocDomain`` section needs to be used as:

.. code-block:: python
   :lineno-start: 55

        RandomLocDomain:
            xmin: 150
            xmax: 250
            ymin: 200
            ymax: 300

This specification will limit the ``x``-coordinate of random wells to be between 150
and 250, and the ``y``-coordinate to be between 200 and 300. Sampling will be
from a uniform distribution on the domain defined by ``xmin`` and ``xmax``
(or ``ymin`` and ``ymax``). Known wells will be placed first; after all known
well coordinates are used, wells will be placed within the random wells domain.
See *ControlFile_ex3.yaml* for an example using random well placement and
*ControlFile_ex4a.yaml* for example using only known well locations.

Note that when using a ``LookupTableReservoir`` component, the wellbore locations
must occur in the domain contained within the files used for the Lookup Table
Reservoir component. If any wellbore locations fall outside of the ``x`` and ``y`` values
covered in that file, the simulation will encounter an error.

Some wellbore components can produce the cumulative masses of |CO2| and brine leaked 
to aquifers (e.g., outputs like **mass_CO2_aquifer1** and **mass_brine_aquifer2**), 
while some cannot. For example, the ``OpenWellbore`` component can produce leakage rates, 
but not leaked masses. In the control file interface, cumulative leaked masses from any 
type of wellbore component will be calculated if ``AccumulateLeakage : True`` is included 
in the wellbore component's section of the *.yaml* control file. These outputs can then 
be saved to ``.csv`` files or used in plots (see *ControlFile_ex48a.yaml*).

Different components can take different entries in a control file. For example, the
``InjectionWell`` entry works with an ``AnalyticalReservoir`` component, but not a
``LookupTableReservoir`` component (for that component, the injection well locations
are set by the reservoir simulations used to create the files used). The entries that
can be provided for each component type are demonstrated in different control file
examples; see chapter :ref:`components_description` to find control file and script
examples for each component.

The last section of the input file is used to specify the plots to be produced.

.. code-block:: python
   :lineno-start: 70

    #-------------------------------------------------
    # Plot setup part of the control file
    #-------------------------------------------------
    Plots:
        CO2_Leakage1:
            TimeSeries: [CO2_aquifer1]
            Subplot:
                NumCols: 2
                Use: True
        CO2_Leakage2:
            TimeSeries: [CO2_aquifer2]
            Subplot:
                NumCols: 2
                Use: True
        Pressure_plot:
            TimeSeries: [pressure]
            Subplot:
                NumCols: 2
                Use: True
                AnalyticalReservoir1_000.pressure: 'Pressure at well #1'
                AnalyticalReservoir1_001.pressure: 'Pressure at well #2'
                AnalyticalReservoir1_002.pressure: 'Pressure at well #3'
                AnalyticalReservoir1_003.pressure: 'Pressure at well #4'
            Title: Reservoir Pressures at Wellbore Locations

Here, three plots are being requested (*CO2_Leakage1*, *CO2_Leakage2*, and
*Pressure_plot*). The first two plots will illustrate the |CO2| leakage to the
shallow aquifer and the thief zone aquifer; the third plot will illustrate the
pressures in the reservoir for the four wellbore locations specified earlier in
the control file. *CO2_Leakage1*, *CO2_Leakage2*, and *Pressure_plot* are the
user-defined names of the three plots to be created; these names will also be used as the
filenames of the figures saved in the output directory. ``TimeSeries`` is a keyword
that instructs the program to plot the observation data as a time series plot. The
values to be plotted in each of the three plots (**CO2_aquifer1**, **CO2_aquifer2**,
and **pressure**) have to be defined in the control file as outputs from one of the
specified component models. Each plot will have a title corresponding to the values
plotted. A user-defined title can be specified with the ``Title`` keyword (as
illustrated for the *Pressure_plot*) in the given plot section. For each aquifer,
the |CO2| leakage rates for all wells will be plotted on the same figure but on
different subplots. If each observation is to be plotted on a separate subplot, the
``Subplot`` keyword with ``Use`` set to ``True`` must be specified, as illustrated
in the example setup. Additionally, the ``NumCols`` keyword (under the ``Subplot``
section) can be used to set the number of subplot columns to use. The number of rows
is controlled by the number of different values (observations) to plot over the number
of columns. Each subplot will be given a default title based on the variable portrayed
in the subplot. The subplot title names can be replaced with the user defined ones by
using the full observation name as a key (e.g., ``AnalyticalReservoir1_000.pressure`` 
for the pressure at the first location) and the desired title as the value under the
``Subplot`` section as shown in the setup of *Pressure_plot*. Note that these titles apply
to specific subplots, while the ``Title`` keyword discussed above applies to the larger figure. 

The example file described here can be found in the *examples/Control_Files* directory with
the filename *ControlFile_ex1a.yaml*. To run this example, open a command prompt in the
*examples/Control_Files* directory, activate the environment created for NRAP-Open-IAM,
and run the command::

    python ../../src/openiam/components/openiam_cf.py --file ControlFile_ex1a.yaml

Note: use \\ on Windows and / on Mac and Linux.

Other example control files can be found in the same directory. These examples can be
run by replacing the file name in the above command with the user specified one. Different 
control file examples are used to demonstrate all of the visualization options available in
NRAP-Open-IAM. For more details regarding the visualization options available, see
section :ref:`cfi_visualization`.

In a control file, an entry situated above the other entries will be read into Python as a
dictionary. For example, the ``Plots`` section shown above will be converted into the
following format in Python:

.. code-block:: python
   :lineno-start: 47

    Plots = {
        'CO2_Leakage1':
            {'TimeSeries': ['CO2_aquifer1'],
            'Subplot':
                 {'NumCols': 2,
                  'Use': True
                  }
            },
         'CO2_Leakage2':
              {'TimeSeries': ['CO2_aquifer2'],
               'Subplot':
                   {'NumCols': 2,
                    'Use': True
                    }
              },
         'Pressure_plot':
             {'TimeSeries': ['pressure'],
              'Subplot':
                 {'NumCols': 2,
                  'Use': True,
                  'AnalyticalReservoir1_000.pressure': 'Pressure at well #1',
                  'AnalyticalReservoir1_001.pressure': 'Pressure at well #2',
                  'AnalyticalReservoir1_002.pressure': 'Pressure at well #3',
                  'AnalyticalReservoir1_003.pressure': 'Pressure at well #4'
                  },
               'Title': 'Reservoir Pressures at Wellbore Locations'
               }
            }

In Python, the '{' and  '}' characters mark the start and end of a dictionary. The entire
control file is read by *openiam_cf.py* as a dictionary called yaml_data. If a user only
wants to run simulations with the CFI, then it is not necessary for the user to understand
this conversion or to know Python. If a user wants to access parts of the CFI through a
script-based approach (e.g., using the NRAP-Open-IAM plot types in a user-created script),
however, then understanding this conversion is important. There are script examples
demonstrating the use of NRAP-Open-IAM plot types developed for the CFI, including the
examples *iam_sys_reservoir_mswell_4aquifers_timeseries.py*,
*iam_sys_reservoir_mswell_stratplot_dipping_strata.py*, and
*iam_sys_reservoir_mswell_futuregen_ttfdplot_dipping_strata.py*.
