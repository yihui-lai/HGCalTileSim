@defgroup PlotUtilsHelper PlotUtilsHelper
@ingroup PlotUtils
@brief Helper function for plot object property computation.

Some of the functions that is needed for generating a nice plot is not directly
provided by the data plot objects. For examples, the TH1D.GetMaximum() returns
the plot range of it's axis, not the actual maximum value of the data inside.
This group of function aims to provide those functions: Common computation
routines that helps make the plot look nicer.
