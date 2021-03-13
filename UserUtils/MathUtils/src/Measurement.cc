/**
 * @file   Measurement.cc
 * @author [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief  Basic functions for Measurement class
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/Measurement.hpp"
#else
#include "UserUtils/MathUtils/Measurement.hpp"
#endif

#include <iostream>

using namespace std;

namespace usr {

/**
 * Default constructor makes an all-zero measurement (may cause issues with
 * calculation... currently un-tested)
 */
Measurement::Measurement()
{
  _central_value = 0;
  _error_up      = 0;
  _error_down    = 0;
}

/**
 * @brief Typical construction involving the assignment of the central value,
 *        the absolute higher uncertainty, and absolute lower uncertainty.
 *
 * Currently the class can only operated with both uncertainties being larger
 * than 0, so if a negative value is passed for the uncertainty, the sign will
 * be automatically flipped.
 */
Measurement::Measurement(
  const double c,
  const double error_up,
  const double error_down ) :
  _central_value( c ),
  _error_up( error_up ),
  _error_down( error_down )
{
  if( _error_up < 0 ){
    cerr << "Warning! Upper error is small than zero! Assuming flipped sign: "
         << c << " " << error_up << " " << error_down
         << endl;
    _error_up = -_error_up;
  }
  if( _error_down < 0 ){
    cerr << "Warning! Lower error is small than zero! Assuming flipped sign: "
         << c << " " << error_up << " " << error_down
         << endl;
    _error_down = -_error_down;
  }
}

/**
 * @brief Assignment operator calls the assignment operator.
 */
Measurement::Measurement( const Measurement& x )
{
  *this = x;
}

/**
 * @brief Nothing to do for the destructor....
 */
Measurement::~Measurement(){}


/**
 * @brief Strict member wise copy constructor.
 */
Measurement&
Measurement::operator=( const Measurement& x )
{
  _central_value = x._central_value;
  _error_up      = x._error_up;
  _error_down    = x._error_down;
  return *this;
}

/*******************************************************************************
*   Arithmetics operators, call functions defined in MeasurementArithmetic.hpp
*******************************************************************************/
Measurement
Measurement::NormParam() const
{
  return Measurement( 1,
    RelUpperError(),
    RelLowerError()
    );
}

/*******************************************************************************
*   Measurement - Measurement arithmetics
*   Call functions defined in MeasurementArithmetic.hpp
*******************************************************************************/
Measurement&
Measurement::operator+=( const Measurement& x )
{
  *this = (*this) + x ;
  return *this;
}

/******************************************************************************/

Measurement&
Measurement::operator*=( const Measurement& x )
{
  *this = (*this) * x ;
  return *this;
}

/******************************************************************************/

Measurement
Measurement::operator+( const Measurement& x ) const
{
  if( x.AbsUpperError() == 0 && x.AbsLowerError() == 0 ){
    return ( *this ) + x.CentralValue();
  } else if( this->AbsUpperError() == 0 && this->AbsLowerError() == 0 ){
    return x + this->CentralValue();
  } else {
    return Sum( *this, x );
  }
}

/******************************************************************************/

Measurement
Measurement::operator*( const Measurement& x ) const
{
  if( x.AbsUpperError() == 0 && x.AbsLowerError() == 0 ){
    return ( *this ) * x.CentralValue();
  } else if( this->AbsUpperError() == 0 && this->AbsLowerError() == 0 ){
    return x + this->CentralValue();
  } else {
    return Prod( *this, x );
  }
}

/*******************************************************************************
*   Measurement - double arithmetics
*   Required for Measurement Arithmetics
*******************************************************************************/
Measurement&
Measurement::operator+=( const double x )
{
  _central_value += x;
  return *this;
}

/*----------------------------------------------------------------------------*/

Measurement&
Measurement::operator-=( const double x )
{
  _central_value -= x;
  return *this;
}

/*----------------------------------------------------------------------------*/

Measurement&
Measurement::operator*=( const double x )
{
  _central_value *= x;
  _error_up      *= fabs( x );
  _error_down    *= fabs( x );
  return *this;
}

/*----------------------------------------------------------------------------*/

Measurement&
Measurement::operator/=( const double x )
{
  _central_value /= x;
  _error_up      /= fabs( x );
  _error_down    /= fabs( x );
  return *this;
}

/*----------------------------------------------------------------------------*/

Measurement
Measurement::operator+( const double x ) const
{
  Measurement ans( *this );
  ans += x;
  return ans;
}
/*----------------------------------------------------------------------------*/

Measurement
Measurement::operator-( const double x ) const
{
  Measurement ans( *this );
  ans -= x;
  return ans;
}

/*----------------------------------------------------------------------------*/

Measurement
Measurement::operator*( const double x ) const
{
  Measurement ans( *this );
  ans *= x;
  return ans;
}

/*----------------------------------------------------------------------------*/

Measurement
Measurement::operator/( const double x ) const
{
  Measurement ans( *this );
  ans /= x;
  return ans;
}

/*----------------------------------------------------------------------------*/

Measurement
operator+( const double y, const Measurement& x )
{
  return x+y;
}

/*----------------------------------------------------------------------------*/

Measurement
operator-( const double y, const Measurement& x )
{
  return Measurement(
    y-x.CentralValue(),
    x.AbsLowerError(),
    x.AbsUpperError()
    );
}

/*----------------------------------------------------------------------------*/

Measurement
operator*( const double y, const Measurement& x )
{
  return x * y;
}

/*----------------------------------------------------------------------------*/

Measurement
operator/( const double y, const Measurement& x )
{
  const double centralValue = y / x._central_value;
  const double err_up       = centralValue * x.RelUpperError();
  const double err_dw       = centralValue * x.RelLowerError();
  return Measurement( centralValue, err_up, err_dw );
}
}/* usr */
