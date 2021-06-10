#include "DataFormats/Common/interface/Wrapper.h"
#include "UserUtils/EDMUtils/interface/Counter.hpp"

namespace ManagerUtil_EDMUtils {
struct dictionaryobjects
{
  usr::Counter               usrct;
  edm::Wrapper<usr::Counter> usrct_wrp;
};
}