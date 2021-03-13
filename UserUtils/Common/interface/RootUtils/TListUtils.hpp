#ifndef USERUTILS_COMMON_ROOTUTILS_TLISTUTILS_HPP
#define USERUTILS_COMMON_ROOTUTILS_TLISTUTILS_HPP

#include "TList.h"

namespace usr {

extern bool MoveObjectToBefore( TList&, const TObject&, const TObject& );

inline bool
MoveObjectToBefore( TList* a, const TObject& b, const TObject& c )
{ return MoveObjectToBefore( *a, b, c ); }

}

#endif
