@class   usr::plt::Pad2DFlat
@ingroup PlotUtilsCore

The creation of this class of to give a uniform way of plotting 2-dimensional
data objects. Since these objects have very distinct methods of representation
from that of 1-dimensional objects (mainly the subjects to be used in in
usr::plt::Pad1D),  a separate class was create to handle the plotting. The
"Flat" in the name of the class indicates that plotting of objects will mainly
focus on projecting information down onto the x-y plane, instead of trying to go
for a 3D projection representation. This representation I find to be easier to
convey and quantify complex diagrams, even if it is not as aesthetically
pleasing for simple data.

A technicality that I feel I need to address is what is the distinction between
"2D" and "1D". Even if scatter plot graphs have data points that contains two
variables (with errors), the data is typically a measurement of one model where
only one of the variables are independent, so I categorized it as "1D". The use
of the "Pad2D" object is where there are two independent variables. There are
cases where one will want to plot a TGraph simple as a collection of points or
as a line on a Pad2D if just for visual guidance. This functionality will be
provided.

The design principle of how this object interacts with ROOT data objects is
similar to the usr::plt::Pad1D object, and will not be repeated here. There are
a couple of distinctions in the plotting behavior that one might need to keep
track of:

- The Legend: Is drawn on the left side of the frame on the outside. We are
  assuming that the entire frame will be filled legends shouldn't over it up.

- The CMS labels, are drawn on the top, outside of the frame.

For the details on the plotting style. See the individual plotting functions.
Support for @ROOT{RooFit} objects have not been implemented yet.
