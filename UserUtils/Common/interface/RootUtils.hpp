#ifndef USERUTILS_COMMON_ROOTUTILS_HPP
#define USERUTILS_COMMON_ROOTUTILS_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/RootUtils/RooArgContainer.hpp"
#include "UserUtils/Common/interface/RootUtils/TListUtils.hpp"
#else
#include "UserUtils/Common/RootUtils/RooArgContainer.hpp"
#include "UserUtils/Common/RootUtils/TListUtils.hpp"
#endif

#endif