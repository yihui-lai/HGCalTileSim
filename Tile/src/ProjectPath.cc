#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/ProjectPath.hh"
#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"
#else
#include "ProjectPath.hh"
#endif

#ifdef CMSSW_GIT_HASH
extern const std::string project_base
  = usr::subpkgpath("HGCalTileSim","Tile" ).string();
#else
extern const std::string project_base = std::string("/data/users/yihuilai/HGCalTileSim") + std::string("/Tile/");
#endif
