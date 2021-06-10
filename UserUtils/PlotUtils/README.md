# `PlotUtils` Package

Helper classes to help maintain a uniform plot style among the plots produced by
analysis code.

## Known issues

- [ ] Simple1DCanvas and Ratio1DCanvas for RooFit cannot be created without a
      second argument (compiler issue? currently testing using CMSSW 9.4.4,
      g++ 6.3.0)
- [ ] Limit support for 1D object:
  - [X] Depreciated: TMultigraph, THStack
  - [ ] TProfile? TEfficient?
- [ ] No support for 2D or higher dimensional objects.
