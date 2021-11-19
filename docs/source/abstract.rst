Abstract
========

The Advanced Circulation Model (ADCIRC) :footcite:p:`Luettich1992` is a Fortran program used for modeling ocean circulation due to tides, surface waves and atmospheric forcings.
However, the input formats and configuration are inflexible and not straight forward for operational implementation, making rapid iteration of model testing, ensemble configuration, and model coupling complicated.
Here, we introduce a flexible abstraction of model inputs and outputs written in Python, called ADCIRCpy :footcite:p:`JaimeCalzadaNOAA2021`, that provides a simpler user interface for automatically generating ADCIRC configuration to a variety of inputs and model scenarios.
This documentation outlines 1. the needs for such an abstraction, 2. the peculiarities and challenges with the ADCIRC model that necessitate custom logic, and 3. methodologies for generalizing user input in such a way as to make generating model configurations consistent, fast, and efficient.

`Python, ADCIRC, configuration generation, circulation modeling, user interface`

.. footbibliography::
