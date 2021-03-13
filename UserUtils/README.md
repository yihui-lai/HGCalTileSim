# UserUtils

This is a package designed to enhance the useability of CMSSW for analysis
purposes and make the related tool sets (ROOT/ROOFIT) more pleasant to use. The
full documentation is hosted [here](https://yimuchen.github.io/UserUtils/).

## Installing

In CMSSW, clone the repository to the `src` directory and compile!

```bash
cd CMSSW_X_Y_Z/src
git clone https://github.com/yimuchen/UserUtils.git
scram b
```

As the plot fixing requires ghostscript, before to have ghostscript installed on
your system. You will also need to export gcc environment variables to allow for
gcc to find the libraries, if this package is not globally installed

```bash
export CFLAGS=$CFLAGS:/your/local/include/path
export LIBRARY_PATH=$LIBRARY_PATH:/your/local/lib/path
```

version requirements: CMSSW environment with gcc version greater than 6 (major
versions 9 and above).

---

As a standalone package, make sure the following is available one your system:

- [cmake](https://cmake.org/) (>3.8)
- [CERN ROOT6](https://root.cern.ch/)
- [boost](https://www.boost.org/) (>1.41)
- [gsl](https://www.gnu.org/software/gsl/)
- [ghostscript](https://www.ghostscript.com/download.html) interpreter interface
- c++ compiler supporting the [c++17 standard](https://en.wikipedia.org/wiki/C%2B%2B17)
  - g++ > 6.0
  - clang > 4

After cloning the repository, build the standalone libraries using `cmake`:

```bash
cd UserUtils
cmake ./
cmake --build ./
```

Currently, this will not install the package as a system library, but place the
required files in directory `UserUtils/lib`. Header files are available in the
`UserUtils/include` directory.

## Code style overview

- The namespace used in this package will always be `usr` if present.
- Header files will always exists in the `interface` directory for the functions
  and classes with the `.hpp` extension. The only exception are EDM plugins
  classes, where a standalone header file will be omitted.
- Source code file will always have the extension of `.cc`. Implementation of
  template functions would placed at the end to their respective header files.
- Class definition, public member functions, and function are in "CamelCase"
  format. Exceptions to this rule include classes that are intended to only be
  used in tandem to STL objects only. In which case they would follow the STL
  convention of all_lowercase_with_underscore.
- No comment is the best comment! I try and keep the control flow as intuitive
  as I can, and keep all comment at the beginning of the function and avoid
  inline comments as much as possible.

## General Known issues/missing functions

See individual packages for package specific issues.

- Testing directories have not been tested locally.

