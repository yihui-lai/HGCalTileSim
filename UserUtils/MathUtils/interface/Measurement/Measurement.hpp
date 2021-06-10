/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Defining the Measurement container class.
 */
#ifndef USERUTILS_MATHUTILS_MEASUREMENT_PARAMETER_HPP
#define USERUTILS_MATHUTILS_MEASUREMENT_PARAMETER_HPP

namespace usr {

/**
 * @brief Class for containing a measurement with asymmetric uncertainties as
 *        well as performing basic arithmetics with quick-and-dirty uncertainty
 *        propagation.
 * @ingroup StatUtils
 */
class Measurement
{
public:
  Measurement();
  Measurement(
    const double central,
    const double error_up,
    const double error_down );
  Measurement( const Measurement& );

  virtual
  ~Measurement ();

  Measurement& operator=( const Measurement& );

  inline operator double() const { return CentralValue(); }

  // Basic access functions
  inline double CentralValue()  const { return _central_value; }
  inline double AbsUpperError() const { return _error_up; }
  inline double AbsLowerError() const { return _error_down; }
  inline double AbsAvgError()   const { return (AbsUpperError()+AbsLowerError() ) / 2; }
  inline double RelUpperError() const { return AbsUpperError()/CentralValue(); }
  inline double RelLowerError() const { return AbsLowerError()/CentralValue(); }
  inline double RelAvgError()   const { return ( RelUpperError()+RelLowerError() ) /2.; }
  inline double UpperValue()    const { return CentralValue() + AbsUpperError();}
  inline double LowerValue()    const { return CentralValue() - AbsLowerError(); }

  // Normalized version
  Measurement NormParam() const;

  // Error arithmetics :: See src/Measurement_Arithmetics.cc
  Measurement&       operator+=( const Measurement& );
  Measurement&       operator*=( const Measurement& );
  Measurement&       operator+=( const double );
  Measurement&       operator-=( const double );
  Measurement&       operator*=( const double );
  Measurement&       operator/=( const double );
  Measurement        operator+( const Measurement& ) const;
  Measurement        operator*( const Measurement& ) const;
  Measurement        operator+( const double )     const;
  Measurement        operator-( const double )     const;
  Measurement        operator*( const double )     const;
  Measurement        operator/( const double )     const;
  friend Measurement operator+( const double, const Measurement& );
  friend Measurement operator-( const double, const Measurement& );
  friend Measurement operator*( const double, const Measurement& );
  friend Measurement operator/( const double, const Measurement& );

private:
  double _central_value;
  double _error_up;
  double _error_down;
};

} /* usr */


#endif/* end of include guard: __MEASUREMENT_HH__ */
