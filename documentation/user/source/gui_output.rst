.. _gui_output:

Output, Plotting, and Analysis in the GUI
=========================================

Output is written to the folder specified in the model definition with the
``Output directory``. If the path to the output directory is not absolute
(i.e., does not containt the drive letter) it is assumed to start from the
NRAP-Open-IAM root folder containing the tool distribution.

For each component model of the system, ``Outputs`` can be specified. When an output
is specified for a forward model the values of that output are written to a file
(*output_name.txt*) in the ``Output Directory``. For a stochastic model (``LHS`` or
``Parstudy``), the parameters used and outputs produced are saved to the files
*parameter_realizations.csv* and *output_realizations.csv*, respectively. A statistical
statistical summary of the input parameters and output observations are also saved
to the files *parameter_stats.csv* and *output_stats.csv*, respectively. Additionally,
the GUI input provided is reformatted into a *.yaml* control file, and this file is
saved in the ``Output Directory`` folder.

After a simulation has run, the post processing section of the GUI can be used
to generate plots of the results and run sensitivity analyses. When the post processing
page is first opened, it asks for a folder specification. This folder is the output folder
containing the results you want to analyze. Navigate to that output folder using the
**Browse** button.

Plotting in the GUI
-------------------
Users can access the post-processing capabilities of GUI by clicking on the
**Post Processing** button on the main page. The **Post Processor** window will
appear.

After a folder containing the results of the simulation is selected, different plotting
options appear. There are several types of plots that can be created depending on the type
of simulation the components used. The simplest plot is a ``Time Series`` plot, where the
output is plotted against time. If the simulation is stochastic (``LHS`` or ``Parstudy``
analysis types), each realization will be plotted as a separate line. When evaluating a stochastic
simulation, the ``Time Series Stats`` and ``Time Series And Stats`` plot types can also be
selected. A ``Time Series Stats`` plot shows basic statistics of the results over time
(mean and median values as well as the four quartiles). A ``Time Series And Stats`` plot
will show these statistics as well as separate lines for each realization.  Use the checkboxes
to select the output types to display in the plot (e.g., **pressure** or **CO2saturation**).
The checkboxes will only reflect the outputs saved by the simulation. The ``Use subplots``
checkbox specifies whether each result should be shown in its own subplot. Finally, the
user can set the plot title and filename with the ``Title`` and ``Filename`` entries,
respectively.

If an ``AtmosphericROM`` component was included in the simulation, map-view plots of the
atmospheric plume can also be generated.

Sensitivity Analysis in the GUI
-------------------------------
If the simulation results are from an ``LHS`` simulation, the ``Processing`` menu
on the ``Post Processing`` window will have options for several types of sensitivity
analysis. Note that while a sensitivity analysis can be run on simulations with a
small number of realizations, the results will most likely be inaccurate. If the sensitivity
coefficients do not sum to one, or if they vary largely through time, the number of
realizations might need to be increased. Generally, 500 to 1000 realizations are needed
for a sensitivity analysis. This number might vary, however, depending on the complexity
of the simulation. Each type of sensitivity analysis will produce plots and/or text file
output in the output directory.

The ``Correlation Coefficients`` option produces a plot matrix of either ``Pearson`` or
``Spearman`` correlation coefficients. Any system model observation can be
excluded from the analysis if needed, although no exclusions need to be made.

The ``Sensitivity Coefficients`` option calculates the sensitivity coefficients
for each selected output to all inputs. Selecting multiple outputs will run
the sensitivity coefficient calculation multiple times. The capture point
is the index for the point in time at which the sensitivity coefficient are
to be calculated. A capture point of 0 will correspond with year 0 of the
simulation. The analysis produces a bar chart.

The ``Multiple Sensitivity Coefficients`` option calculates the impact of input parameters
on multiple outputs. Multiple outputs should be selected here. The capture point
is the index for point in time at which the sensitivity coefficient are to be calculated.
The analysis will produce a bar chart.

The ``Time Series Sensitivity`` option will produce a line graph illustrating how the impact
from input parameters changes over time with respect to an output value. Selecting multiple
output values will run the analysis multiple times. The capture point determines the time
at which the sensitivity coefficients are compared and then ordered based on the comparison.
