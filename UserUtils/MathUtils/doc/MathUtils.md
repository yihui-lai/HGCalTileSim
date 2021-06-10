@defgroup MathUtils MathUtils
@brief Collection of common mathematical calculation routines used in HEP
       analysis.

This is a collection of functions and classes used for the mathematical routines
commonly used in the HEP analysis that are not immediately available in existing
libraries, the main ones being the calculation of uncertainties of measurements
with arithmetics, and the Kolmogorov-Smirnov test weighted data sets. This
libraries uses custom written GSL functions for the bulk of the calculation, but
will use existing ROOT functions if they are available.
