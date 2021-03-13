/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   New formats to allow for measurement/ostream interaction
 */
#ifndef USERUTILS_MATHUTILS_MEASUREMENT_FORMAT_HPP
#define USERUTILS_MATHUTILS_MEASUREMENT_FORMAT_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/BoostUtils/PTreeUtils.hpp"
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"
#else
#include "UserUtils/Common/BoostUtils/PTreeUtils.hpp"
#include "UserUtils/Common/Format.hpp"
#include "UserUtils/MathUtils/Measurement/Measurement.hpp"
#endif

#include <string>

#include "RooRealVar.h"

namespace usr {

namespace fmt {

/**
 * @brief decimal representation of a measurement with uncertainties.
 * @ingroup StatUtils
 */
class decimal : public base::format
{
public:
  decimal(
    const double input,
    const int    precision = -base::precision_default );
  decimal(
    const Measurement& input,
    const int          precision = -base::precision_default );
  decimal(
    const RooRealVar& input,
    const int         precision = -base::precision_default );
  decimal(
    const RooRealVar* input,
    const int         precision = -base::precision_default );

  virtual
  ~decimal (){}

  inline decimal& precision( const int x );// redefining

  virtual std::string str() const;// Explicit function call

private:
  double _central;
  double _upper;
  double _lower;

  void SetPrecision();
};

/**
 * @brief Scientific representation of a measurement with uncertainties.
 * @ingroup StatUtils
 */
class scientific : public base::format
{
public:
  scientific(
    const double   input,
    const unsigned precision );
  scientific(
    const Measurement& input,
    const int          precision = -1 ); // default to auto for
  scientific(
    const RooRealVar& input,
    const int         precision = -1 );
  scientific(
    const RooRealVar* input,
    const int         precision = -1 );
  virtual
  ~scientific  (){}

  inline scientific& precision( const int x );// redefining

  std::string str() const;// Explicit function call

private:
  double _central;
  double _upper;
  double _lower;
  unsigned _exp;

  // Private helper functions
  void SetExponent();
  void SetPrecision();
};

}/* fmt */


/*-----------------------------------------------------------------------------
 *  Additional IO formatting options
   --------------------------------------------------------------------------*/
template<>
Measurement GetSingle<Measurement>(
  const pt::ptree&,
  const std::string&
  );

}/* usr */

#endif/* end of include guard: USERUTILS_MATHUTILS_PARAMETERFORMAT_HPP */
