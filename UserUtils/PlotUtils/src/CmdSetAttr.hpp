#ifndef USERUTILS_PLOTUTILS_SRC_CMDSETATTR_HPP
#define USERUTILS_PLOTUTILS_SRC_CMDSETATTR_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/RootUtils/RooArgContainer.hpp"
#else
#include "UserUtils/Common/RootUtils/RooArgContainer.hpp"
#endif

#include "TAttText.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttMarker.h"

extern void SetTextAttr( TAttText& obj,
  const usr::RooArgContainer& args  );

extern void SetLineAttr( TAttLine& obj,
  const usr::RooArgContainer& args  );

extern void SetFillAttr( TAttFill& obj,
  const usr::RooArgContainer& args  );

extern void SetMarkAttr( TAttMarker& obj,
  const usr::RooArgContainer& args  );

#endif
