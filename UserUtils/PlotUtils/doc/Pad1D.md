@class   usr::plt::Pad1D
@ingroup PlotUtilsCore

This class aims for a uniform plotting interface for both @ROOT and @ROOFIT
objects, assisting with plotting objects according to CMS recommended plotting
style, removing the use of 'magic strings' needed for object plotting, reducing
the verbosity of the code needed for generating publication ready plots, and
improving the readability of the plotting related codes in general.

# Plotting ROOT and RooFit objects under the same interface

@ROOT and @ROOFIT takes two distinct style of plotting: In @ROOT, data
containers and plotting objects are one-and-the-same (i.e. a @ROOT{TH1D}
instance is used for data interaction via `Fill` functions, and stylising the
plotting output via the `SetLineColor()` and similar functions); In @ROOFIT, the
data containers (@ROOT{RooDataSet}s, and @ROOT{RooAbsPdf}) are used to generated
plotting classes (@ROOT{TGraph}s) via the `plotOn` functions, with the
assistance of a @ROOT{RooPlot} object, where the generated objects are destroyed
when the @ROOT{RooPlot} object goes out of scope.

@ROOFIT has a very good reason for doing such a thing: the process of generating
a simple plotting object from a complicated, multi-variable entity is
non-trivial, and in the case of serious numerical analysis, the same data entity
might be required to be plotted multiple times with different extraction
methods. Defining the plotting styles in the objects could very likely cause
unwanted behaviors in the resulting plots. However, @ROOT objects are relatively
simple, with very few ways of presenting the stored data, so not using an
additional generated object to plot and object 'makes sense'. The final design
decision was to keep the plotting of @ROOT objects "ROOT"-like; and @ROOFIT
objects "RooFit"-like:

- @ROOT objects are drawn on the Pad instance, with the objects ownership
  untouched. A wrapper for the `PadBase::PlotObj()` functions would be used for
  additional style handling.
- @ROOFIT object would have their representative @ROOT{TGraph} objects generated
  by an internal @ROOT{RooPlot} object, then the @ROOT{TGraph} will have their
  plotting options passed onto as to the ROOT object plotting functions. This is
  to avoid repeated calls to the `RooPlot::Draw()` functions, as this function
  basically loops overall the internal objects and passes it to the working
  @ROOT{TPad}. The @ROOFIT plotting functions would also return a reference to
  the generated @ROOT{TGraph} object to allow for additional stylising.

The plotting functions will essentially be a improved @ROOFIT plotting interface
using the @ROOT{RooCmdArg} class to change the behavior of the plotting:

```cpp
pad.Plot<ObjType>( PlotObj, Cmd1, Cmd2, Cmd3,... );
```

For the case of plotting @ROOFIT objects, all RooFit defined @ROOT{RooCmdArg}s
are acceptable, with some custom commands made designed additional enhancements.

# Object plotting behavior enhancements

A few addition plotting behavior for objects are present:

- **PlotType** allows for the specification of how the data should be presented:
  whether as scatter plot (data point with error bars), a function line with
  uncertainty band or standard solid histograms or histogram error boxes.
  Additional parsing for setting the styling of with CMS plotting conventions
  (ex. no x error bars should be present for equal-bin-width histograms...etc)
  are also present. More details would be given in the Pad1D::PlotHist() and
  Pad1D::PlotGraph() functions.

- **TrackY** tells the pad whether for not to adjust the Y axis range
  depending on input object. If it is not specified, it is automatically
  determined as followed:
  - If the object is a histogram, then only the maximum value is tracked for
    axis range adjustment (the default 'minimum value' is 0.3).
  - If the object is a @ROOT{TGraph}, and no histogram is already present in the
    pad, then both the maximum and minimum values are tracked. If any histogram
    is already present however, the new @ROOT{TGraph} would not change the axis
    range.

- **EntryText** tells the object to be added to the legend on the pad. What is
  to be presented in the legend would be automatically determined by the
  PlotType passed to the plot function.

- **VisualizeError** is now extended to include TF1, which takes a TFitResult to
  interpret for the error band should be drawn.

- **PlotUnder**. In ROOT, the first plotted objects are always placed at the back
  of the canvas. Sometimes this is not desireable, since the a object that is
  plotted later might want to be placed in the background, while required to be
  calculated only after a certain object has been generated (an example that
  comes to mind is the plotting of the RooFit objects, where the plotting of
  @ROOT{RooAbsPdf}s depends on the existence of a @ROOT{RooAbsData} object,
  unless explicitly specified). This option will be force a later plotted object
  to be behind a certain target.

# Additional Pad plotting behavior enhancements

- Standard CMS and luminosity labels text would be added to the pad via the
  Pad1D::DrawCMSLabel() and Pad1D::DrawLuminosity() functions. The CMS label
  would be placed on the top-left side inside the axis frame, while the
  luminosity label would be placed on the top-right side outside the axis frame.
- If an EntryText is specified at any point in the pad's lifetime, A
  @ROOT{TLegend} would be spawn on the the pad inside the axis frame. The size of
  the frame would be automatically adjusted, based on the entry text lengths and
  the number of entries specified. The position of the legend can be any of the 9
  positions in the left--center--right, top--center--bottom system. By default,
  the legend is placed on the top-right, exploiting the fact that HEP 1D plots
  typically exhibit long trailing tails on the right.

  An issue with ROOT object drawing is that the first object to be
  drawn has a smaller z index (placed towards the back), when usually the
  plotted towards the foreground are more 'important' (should be higher on
  the legend entries). So we will be using a stack to temporary store
  the legend entry and flushing it's content when a save function is called.
  When the flushing is called, additional method would also be used to
  determine the size of the legend, and the size of the legend would be
  adjusted accordingly.
