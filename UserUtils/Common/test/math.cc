#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/Common/interface/Maths.hpp"
#else
#include "UserUtils/Common/Format.hpp"
#include "UserUtils/Common/Maths.hpp"
#endif

#include <iostream>

int main()
{
  uint32_t a = usr::HashValue32( 0.123     );
  uint32_t b = usr::HashValue32( 0.1230001 );
  uint32_t c = usr::HashValue32( 0.12300000000000001 );
  uint32_t d = usr::HashValue32( 100.00 );

  std::cout << a << " " << b << " " << c << std::endl;
  std::cout << usr::Hash32Join( a, b ) << " "
            << usr::Hash32Join( b, a ) << std::endl;

  std::cout << usr::Hash32Join( usr::Hash32Join( a, b ), a ) << " "
            << usr::Hash32Join( usr::Hash32Join( a, a ), b )  << std::endl;
  return 0;
}
