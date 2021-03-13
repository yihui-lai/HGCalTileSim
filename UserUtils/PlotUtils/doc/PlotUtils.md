@defgroup PlotUtils PlotUtils
@brief Better plotting functions and classes.

This subpackage is a collections of classes and functions intended to enhance
the plotting experience for analysis level plots by augmenting the existing
@ROOT classes.

Due to the complexity of the issues with root. This package is further divided
into the following subpackages:

- [PlotUtilsUser](@ref PlotUtilsUser)
  User level classes to allow for presentation ready plots with minimal user
  code.

- [PlotUtilsArgument](@ref PlotUtilsArgument)
  List of arguments that can be used for plotting.

- [PlotUtilsHelper](@ref PlotUtilsHelper)
  Functions to help with plot object property computations.

- [PlotUtilsCore](@ref PlotUtilsCore):
  Are the key wrapper classes that *restricts* the usage of ROOT plotting
  objects, wrapping the plotting APIs of ROOT into more intuitive methods,
  which hopefully improves code readability. Also includes classes for
  conventional styling conventions used by the CMS collaboration.

- [PlotUtilsConstants](@ref PlotUtilsConstants)
  Translation of the many magic numbers used in ROOTs plotting API to human
  readable constants and functions.
