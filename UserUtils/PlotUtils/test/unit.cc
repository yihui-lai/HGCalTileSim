/**
 * @file    unit_test.cc
 * @brief   Test for unit conversion
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Constants.hpp"
#include "UserUtils/Common/interface/STLUtils.hpp"
#else
#include "UserUtils/PlotUtils/Constants.hpp"
#include "UserUtils/Common/STLUtils.hpp"
#endif

#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
  cout << usr::separator() << endl
       << "Length conversion and constants" << endl;
  {
    cout << usr::plt::len::a4paperwidth() << endl
         << usr::plt::len::a4textheight_default() << endl;

    cout << usr::plt::len::a4textwidth_default() << endl
         << usr::plt::len::a4textheight_default() << endl;
  }

  return 0;
}
