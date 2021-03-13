/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   One file include for all STL utility headers
 */
#ifndef USERUTILS_COMMON_STLUTILS_HPP
#define USERUTILS_COMMON_STLUTILS_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"
#else
#include "UserUtils/Common/STLUtils/OStreamUtils.hpp"
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/STLUtils/Filesystem.hpp"
#endif

#endif /* end of include guard: USERUTILS_COMMON_STLUTILS_HPP */
