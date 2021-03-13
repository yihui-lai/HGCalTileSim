#!/bin/bash
# Helper script for generating the the reflection CPP files for python binding

HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BINDING_PKG="MathUtils PlotUtils"

for PKG in $BINDING_PKG; do
  DIR_NAME=$HOME_DIR/${PKG}/plugins/
  HEADER=${PKG}.hpp
  SELECTION_NAME=${PKG}.xml
  OUTPUT_CC=${PKG}_rflx.cc
  cd   $DIR_NAME
  genreflex -I$HOME_DIR/include \
            -s $SELECTION_NAME  \
            -o $OUTPUT_CC       \
            $HEADER
done
