.. _cfi_analysis:

Analysis Options in the CFI
===========================

NRAP-Open-IAM uses the Model Analysis ToolKit (MATK) :cite:`MATK` for
the basis of its probabilistic framework. More information about MATK can be found here:
`http://dharp.github.io/matk/ <http://dharp.github.io/matk/>`_. The MATK code repository can
be found here: `https://github.com/dharp/matk <https://github.com/dharp/matk>`_.

Parameter input and output interactions can be explored using the *Analysis* section
of a *.yaml* control file. Correlation coefficients can be calculated using the
``CorrelationCoeff`` keyword. Parameter sensitivity coefficients for any output
simulation value can be calculated using a Random-Balanced-Design Fourier
Amplitude Sensitivity Test (RBD-Fast) technique. The control file keywords
``SensitivityCoeff``, ``MultiSensitivities``, ``TimeSeriesSensitivity`` can be
used to access different sensitivity coefficient outputs. See control file examples 
*ControlFile_ex8a.yaml* to *ControlFile_ex8d.yaml* for more details regarding the
use of the analysis section.

The Sensitivity Analysis is done with the ``SALib`` package :cite:`Herman2017`.
For more information on the RBD-Fast technique see :cite:`TISSOT2012205` and
:cite:`PLISCHKE2010354`. While not accessible through the control files, a
scripting interface is provided to a Sobol sensitivity analysis :cite:`SOBOL20093009`.
