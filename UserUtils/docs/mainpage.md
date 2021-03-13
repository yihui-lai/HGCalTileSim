## What is UserUtils?

UserUtils started out as random packages I wrote from when I was completing my
masters thesis. CMSSW is great for organising analysis code, but has a very
sparse documentation for recipes, when taking into CMSSW framework related
recipes, analysis related recipes, STL/boost recipes, ROOT/RooFit recipes are
all required for analysis work. Later, after I finished my masters theisis, I
clean up and organized the many functions and class I wrote into this package:
a collection of general purpose classes and functions aimed to enhance the
programming experience for HEP analysis. Hopefully, even if you don't directly
use this package, this could serve as a starting point to get you started with
HEP analysis codes.

## Contents

As I was programming in the CMSSW framework, and using the data format designed
by the CMS computing team, not all functions would work outside the CMSSW
framework. I have, however, managed to split the packages that doesn't
explicitly depend on the CMSSW framework into standalone packages that could be
compiled independently. As long as I still work in CMS however, those packages
would be maintained.

Installation instructions and requirements could be found on the source code
[repository](https://github.com/yimuchen/UserUtils)

### Standalone packages

* **[Common](@ref Common)**:
  A group of STL and boost related recipes such as full-range vector operations,
  better floating point `ostream` formatting, simplified interface for
  `boost::property_tree` interactions... etc.

* **[MathUtils](@ref MathUtils)**:
  Mathematical functions that I used in my analysis that I consider "general
  purpose" but is not directly implemented in any of the libraries
  (GSL/ROOT/RooFit). This includes a container for a measurement with
  uncertainties along with quick-and-dirty uncertainty propagation with
  arithmetic operations, and goodness-of-fit statistics with RooFit objects.

* **[PlotUtils](@ref PlotUtils)**:
  Overhaul of the plotting interface of ROOT, with a better wrapper class for
  Canvases to allow for a more intuitive interface for analysis level plotting.

### CMSSW packages

* **[EDMUtils](@ref EDMUtils)**:
  Various EDM utilities including classes to allow for caching of EDM file
  contents, alias classes for reducing plugin verbosity, scripts for dumping
  specific EDM file contents, and standalone plugins for common information
  caching.

* **[PhysUtils](@ref PhysUtils)**:
  Recipes for interacting with the physics objects in CMSSW, including
  trigger matching with glob notation, and generation truth topology crawling,
  and various additional variable calculations.

---

## Contributing

The source code is hosted on [Github](https://github.com/yimuchen/UserUtils),
so please feel free to contribute!
