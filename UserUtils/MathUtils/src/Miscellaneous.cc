/**
 * @file    Miscellaneous.cc
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief  Implementing miscellaneous mathematical functions
 */

#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/Miscellaneous.hpp"
#else
#include "UserUtils/MathUtils/Miscellaneous.hpp"
#endif
#include "TMath.h"

#include <cmath>

namespace usr {

/**
 * @brief returning the intersect points of two segments if exists.
 *
 * Given four points \f$A(p1x,p1y)\f$, \f$B(p2x,p2y)\f$, \f$C(p3x,p3y)\f$,
 * and \f$D(p4x,p4y)\f$. The function checks if segments AB and CD intersect.
 * If yes, store results in intx, and inty return true; return false otherwise.
 */
bool
Intersect(
  const double p1x, const double p1y,
  const double p2x, const double p2y,
  const double p3x, const double p3y,
  const double p4x, const double p4y,
  double& intx, double& inty
  )
{
  const double s1x = p2x - p1x;
  const double s1y = p2y - p1y;
  const double s2x = p4x - p3x;
  const double s2y = p4y - p3y;

  const double s = ( -s1y * ( p1x - p3x ) + s1x * ( p1y - p3y ) ) / ( -s2x*s1y + s1x*s2y );
  const double t = ( s2x * ( p1y - p3y ) - s2y * ( p1x - p3x ) ) / ( -s2x*s1y + s1x*s2y );

  if( s >= 0 && s <= 1 && t >= 0 && t <= 1 ){
    intx = p1x + ( t * s1x );
    inty = p1y + ( t * s1y );
    return true;
  } else {
    return false;
  }
}

/**
 * @brief Performing a Cholesky decomposing of a covariance matrix from fit
 *
 * In the case that the input matrix is not positive definite, this function will
 * a diagonal matrix  that is just the square root of the diagonal components of
 * the original covariance matrix. The function is essentially the
 * [TDecompChol::Decompose](https://root.cern.ch/doc/master/classTDecompChol.html#af14df10a3c766330cb93063161ecedc0) method.
 */
extern TMatrixD
DecompCorvariance( const TMatrixDSym& m )
{
  const unsigned n = m.GetNrows();
  TMatrixD ans     = m;
  double* ans_ptr  = ans.GetMatrixArray();

  for( unsigned icol = 0; icol < n; icol++ ){
    const unsigned rowOff = icol*n;

    double ujj = ans_ptr[rowOff+icol];

    for( unsigned irow = 0; irow < icol; irow++ ){
      const unsigned pos_ij = irow*n+icol;
      ujj -= ans_ptr[pos_ij]*ans_ptr[pos_ij];
    }

    if( ujj <= 0 ){// Special returning case for any non-definite positive case
      for( unsigned i = 0; i < n; ++i ){
        const unsigned offset = i * n;

        for( unsigned j = 0; j < n; ++j ){
          ans_ptr[i+offset] = i != j ? 0 :
                              TMath::Sqrt( m[i][j] );
        }
      }

      return ans;
    }
    ujj                  = TMath::Sqrt( ujj );
    ans_ptr[rowOff+icol] = ujj;

    if( icol < n-1 ){
      for( unsigned j = icol+1; j < n; j++ ){
        for( unsigned i = 0; i < icol; i++ ){
          const unsigned rowOff2 = i*n;
          ans_ptr[rowOff+j] -= ans_ptr[rowOff2+j] * ans_ptr[rowOff2+icol];
        }
      }

      for( unsigned j = icol+1; j < n; j++ ){
        ans_ptr[rowOff+j] /= ujj;
      }
    }
  }

  for( unsigned irow = 0; irow < n; irow++ ){
    const unsigned rowOff = irow*n;

    for( unsigned icol = 0; icol < irow; icol++ ){
      ans_ptr[rowOff+icol] = 0.;
    }
  }

  return ans;
}

}/* usr */
