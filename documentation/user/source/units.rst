Units
=====

Data passed between models need to have consistent units. Here
is a list of the units used by NRAP-Open-IAM.

* Pressure is assumed to be in units of Pascals (|Pa|).
* Time is assumed to be in days.
* Distance, width, length, height, depth are assumed to be in units of meters (|m|).
* Flow rates are assumed to be in units of kilograms per second (|kg/s|).
* Mass is assume to be in units of kilograms (|kg|).
* Viscosities are assumed to be in units of Pascal seconds (|Pa*s|).
* Permeability is assumed to be in units of meters squared (|m^2|).

Note that these rules apply to the exchange of data between component models in
NRAP-Open-IAM during a simulation. There are input options in the CFI and GUI that
take values in different units. For example, the ``AoR`` plot and ``AoR`` workflow
have the input ``CriticalPressureMPa``, which takes a pressure value in |MPa|. While
NRAP-Open-IAM simulations use time in days, many temporal inputs provided through
the CFI or GUI take time in years (before starting the simulation, the conversion to
days is handled automatically). Any such distinctions are clearly noted in this document.
