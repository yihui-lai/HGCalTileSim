/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Implementing measurement ostream interaction.
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement.hpp"
#else
#include "UserUtils/Common/Maths.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/Measurement.hpp"
#endif

#include <boost/algorithm/string.hpp>
#include <regex>
#include <string>

namespace usr
{

namespace fmt {

/**
 * @brief Re-implementing a double input for common fmt interface
 *
 * This construction means that the measurement without uncertainty would
 * output the same values as if using a double.
 */
decimal::decimal( const double input, const int p ) :
  _central( input ),
  _upper( 0 ),
  _lower( 0 )
{
  precision( p );
}

/**
 * @brief normal construction for typical measurements.
 */
decimal::decimal( const Measurement& input, const int p ) :
  _central( input.CentralValue() ),
  _upper( input.AbsUpperError() ),
  _lower( input.AbsLowerError() )
{
  if( p < 0 ){
    SetPrecision();
  } else {
    precision( p );
  }
}

/**
 * @brief Construct using a RooRealVar status.
 */
decimal::decimal( const RooRealVar& input, const int p ) :
  _central( input.getVal() ),
  _upper( fabs( input.getErrorHi() ) ),
  _lower( fabs( input.getErrorLo() ) )
{
  if( p < 0 ){
    SetPrecision();
  } else {
    precision( p );
  }
}

/**
 * @brief Construct using a RooRealVar status.
 */
decimal::decimal( const RooRealVar* input, const int p ) :
  _central( input->getVal() ),
  _upper( fabs( input->getErrorHi() ) ),
  _lower( fabs( input->getErrorLo() ) )
{
  if( p < 0 ){
    SetPrecision();
  } else {
    precision( p );
  }
}

/**
 * @brief main operation for creating a latex string.
 *
 * @details Note that this function essentially calls the base::decimal
 * methods for creating the string representation of a double three times,
 * meaning that if the central value and uncertainties are of wildly different
 * orders of magnitude, the output could look very weird.
 */
std::string
decimal::str() const
{
  if( _upper == _lower && _upper == 0 ){
    return base::decimal( _central ).dupsetting( *this ).str();
  } else {
    const std::string cen = base::decimal( _central ).dupsetting( *this ).str();
    const std::string up  = base::decimal( _upper ).dupsetting( *this ).str();
    const std::string lo  = base::decimal( _lower ).dupsetting( *this ).str();

    if( up == lo ){// string-wise comparison!
      return usr::fstr( "%s\\pm%s", cen, up );
    } else {
      return usr::fstr( "%s^{+%s}_{-%s}", cen, up, lo );
    }
  }
}

/**
 * @brief Base precision overloaded to allow for autoamtic precision setting
 *        with negative settings.
 */
decimal&
decimal::precision( const int p )
{
  if( p < 0 && ( _upper != 0 || _lower != 0 ) ){
    SetPrecision();
  } else {
    base::format::precision( abs( p ) );
  }
  return *this;
}

/**
 * @brief automatic precision setting for the decimal representation.
 * @details The precision is set such that the larger uncertain would display
 *          at least two significant digits (0 if larger uncertainty is greater
 *          than 1).
 */
void
decimal::SetPrecision()
{
  const double op_unc = std::max( _upper, _lower );
  const int exp       = GetExponent( op_unc );
  if( exp > 0 ){
    _precision = 0;
  } else {
    _precision = -exp+1;
  }
}

/*-----------------------------------------------------------------------------
 *  Scientific notation implementation functions
   --------------------------------------------------------------------------*/


/**
 * @brief allowing double input to avoid multiple interfaces for double input.
 *
 * This construction means that the measurement without uncertainty would
 * output the same values as if using a double.
 */
scientific::scientific( const double input, const unsigned p ) :
  _central( input ),
  _upper( 0 ),
  _lower( 0 ),
  _exp( 0 )
{
  precision( p );
  SetExponent();
}

/**
 * @brief simple construction with a measurement with uncertainties.
 */
scientific::scientific( const Measurement& input, const int p ) :
  _central( input.CentralValue() ),
  _upper( input.AbsUpperError() ),
  _lower( input.AbsLowerError() ),
  _exp( 0 )
{
  SetExponent();
  if( p < 0 ){
    SetPrecision();
  } else {
    precision( p );
  }
}

/**
 * @brief Construct using a RooRealVar instance
 */
scientific::scientific( const RooRealVar& input, const int p ) :
  _central( input.getVal() ),
  _upper( fabs( input.getErrorHi() ) ),
  _lower( fabs( input.getErrorLo() ) ),
  _exp( 0 )
{
  SetExponent();
  if( p < 0 ){
    SetPrecision();
  } else {
    precision( p );
  }
}

/**
 * @brief Construct using a RooRealVar instance
 */
scientific::scientific( const RooRealVar* input, const int p ) :
  _central( input->getVal() ),
  _upper( fabs( input->getErrorHi() ) ),
  _lower( fabs( input->getErrorLo() ) ),
  _exp( 0 )
{
  SetExponent();
  if( p < 0 ){
    SetPrecision();
  } else {
    precision( p );
  }
}

/**
 * @brief implementing the virtual function.
 *
 * The function attempts to use the simplest form to represent the uncertainty,
 * using as little latex symbols as possible.
 */
  std::string
  scientific::str() const
{
  const std::string cen  = base::decimal( _central ).dupsetting( *this ).str();
  const std::string up   = base::decimal( _upper ).dupsetting( *this ).str();
  const std::string lo   = decimal( _lower ).dupsetting( *this ).str();
  const std::string base =
    ( up == lo && _upper == 0 ) ? cen :
    ( up == lo )                ? usr::fstr( "%s\\pm%s", cen, up ) :
    usr::fstr( "%s^{+%s}_{-%s}", cen, up, lo );

  const std::string ans =
    ( _exp == 0 )               ? base :
    ( up == lo && _upper != 0 ) ? usr::fstr( "(%s)\\times10^{%d}", base, _exp ) :
    usr::fstr( "%s\\times10^{%d}", base, _exp );

  return ans;
}

/**
 * @brief overloading base implementation to allow for precision autodetection
 *        with a negative input.
 */
scientific&
scientific::precision( const int i )
{
  if( i < 0 && ( _upper != 0 || _lower != 0 ) ){
    SetPrecision();
  } else {
    base::format::precision( abs( i ) );
  }
  return *this;
}

/**
 * @brief reducing/magnifying the central value and uncertainties by a common
 *        exponent value.
 *
 * If the central value is none-zero, then extracting the exponent of the
 * central value such that \f$1 < |\mathrm{central}_\mathrm{man}| < 10\f$;
 * If the central value is zero, then the exponent is determined by the
 * larger of the uncertainties, if it is non-zero.
 * If everything is zero, then the exponent is not set.
 */
void
scientific::SetExponent()
{
  if( _central != 0 ){
    _exp      = GetExponent( _central );
    _central /= IntPower( 10, _exp );
    _upper   /= IntPower( 10, _exp );
    _lower   /= IntPower( 10, _exp );
  } else if( _upper != 0 || _lower != 0  ){
    _exp      = std::max( GetExponent( _upper ), GetExponent( _lower ) );
    _central /= IntPower( 10, _exp );
    _upper   /= IntPower( 10, _exp );
    _lower   /= IntPower( 10, _exp );
  }
}

/**
 * @brief automatically determining the precision to use.
 *
 * If the central value is smaller than the larger uncertainty, then at most
 * the precision is 1 (if the uncertainty after the exponent has been factored
 * out smaller than 10).
 * if the central value is larger than the larger uncertainty, the precision
 * is set such that two significant digits of the larger uncertainty is
 * displayed.
 */
void
scientific::SetPrecision()
{
  double op_unc = std::max( _upper, _lower );
  if( op_unc > _central ){
    if( op_unc < 10 ){
      _precision = 1;
    } else {
      _precision = 0;
    }
  } else {
    _precision = 1;

    while( op_unc < 1 ){
      _precision++;
      op_unc *= 10;
    }
  }
}

}/* fmt */


/**
 * @brief template specialization for a boost property tree helper function.
 *
 * reading a list of doubles with 1-3 parameters as a measurement class.
 * * If only an empty list exists at the address, then the unit measurement
 *   (1+-0), * is returned.
 * * If only the list is only a single double long, a measurement with no
 *   uncertainty is returned.
 * * If only two double exist, then a measurement with symmetric uncertainty is
 *   returned. (x[0] +- x[1])
 * * If more than three doubles exists, then the first three doubles will be
 *   used to construct a measurement with asymmetric uncertainty:
 *   (x[0] +x[1] -x[2] )
 */
template<>
Measurement
GetSingle<Measurement>(
  const pt::ptree&   tree,
  const std::string& query )
{
  const std::vector<double> input = GetList<double>( tree, query );
  if( input.size() == 0 ){
    return Measurement( 1, 0, 0 );// Returning unit for empty input
  } else if( input.size() == 1 ){
    return Measurement( input.at( 0 ), 0, 0 );
  } else if( input.size() == 2 ){
    return Measurement( input.at( 0 ), input.at( 1 ), input.at( 1 ) );
  } else {
    return Measurement( input[0], input[1], input[2] );
  }
}

}/* usr */
