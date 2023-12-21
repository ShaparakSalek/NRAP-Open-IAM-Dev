.. _equations:

Equations
=========

In this section, we present equations related to geologic carbon storage. While all component
models have their own governing equations, information about each component can be found in chapter
:ref:`components_description`. This section is only meant to address equations that are frequently
referenced in this document.

Critical Pressure
-----------------

Critical pressure can be used in ``AoR`` plots as well as the leakage calculations of the ``OpenWellbore``
component. If given a critical pressure, an ``AoR`` plot evaluating reservoir pressures will highlight any
points where the pressures exceed the critical pressure. In this case, the ``AoR`` plot is only considering
the output of a reservoir component (not a wellbore or aquifer component).

If the critical pressure (|Pcrit|) is not given as a specific value, it is calculated as :cite:`USEPA2013`:

    |Pcrit| = (|rho_w| |times| g |times| |d_aq|) + (|rho_br| |times| g |times| (|d_res| - |d_aq|)),

where |rho_w| and |rho_br| are the densities of water (1000 |kg/m^3|) and brine, respectively, 
g is gravitational acceleration (9.81 |m/s^2|), |d_aq| is the depth to the bottom of the aquifer impacted 
by leakage (|m|), and |d_res| is the depth to the top of the reservoir (|m|). Higher brine densities generally
produce lower leakage rates; it is more difficult to lift brine from the reservoir to an aquifer when the
brine has a higher density. In an ``AoR`` plot evaluating reservoir pressure, |d_aq| reflects the aquifer being 
considered in the ``AoR`` analysis. For an ``OpenWellbore`` component using a critical pressure, |d_aq|
should reflect the aquifer receiving leakage from the wellbore (i.e., |d_aq| is equal to the **wellTop**
parameter) and |rho_br| is set by the **brineDensity** parameter.

If an ``AoR`` plot is made for a simulation using a wellbore component that has a brine density parameter
(e.g., **brineDensity**  parameter of the ``OpenWellbore`` and ``MultisegmentedWellbore`` components),
then that parameter will be used for the |rho_br| value. If the wellbore component does not have a brine
density parameter (e.g., ``CementedWellbore``), then a brine density can be specified with the ``BrineDensity``
entry for the ``AoR`` plot type. This entry can also be used with the ``AoR`` workflow. Otherwise, a default
brine density value of 1045 |kg/m^3| will be used.

