/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Classes for better floating point string/ostream outputs
 */
#ifndef USERUTILS_COMMON_FORMAT_HPP
#define USERUTILS_COMMON_FORMAT_HPP

#include <ostream>
#include <string>

namespace usr {

namespace fmt {

namespace base {

/**
 * @addtogroup ostreamutils
 * @{
 */

extern int precision_default;
extern unsigned spacesep_default;
extern std::string spacestr_default;
extern const unsigned max_precision;

/**
 * @brief Base class for defining float formatting operations.
 * @details Pretty much an all settings class. with all settings function
 *          returning 'this reference' to allow for setting chaining.
 */
class format
{
public:
  virtual ~format(){};

  /**
   * @brief number of digits to display after the decimal point.
   */
  inline format&
  precision( int x ){ _precision = x; return *this; }

  /**
   * @brief string to add between digit seperations.
   */
  inline format&
  spacestr( const std::string& x ){ _spacestr = x; return *this; }

  /**
   * @brief number of digits between seperator string.
   */
  inline format&
  spacesep( unsigned x ){ _spacesep = x; return *this; }


  format& operator=( const format& ) = default;// Member-wise copy

  /**
   * @brief duplicating all settings to reference format.
   */
  inline format&
  dupsetting( const format& x ){ *this = x; return *this;}


  /**
   * @brief Pure virtual function for determining how an type should be
   * represented as a string.
   */
  virtual std::string str() const = 0;

  /**
   * @brief conversion operator to string.
   */
  inline operator std::string(){
    return str();
  }

protected:
  int _precision;
  std::string _spacestr;
  unsigned _spacesep;

  /**
   * @brief constructor can only be called by childern classes.
   * @details intentionally hidden to stop user usage.
   */
  format() : _precision( precision_default ),
    _spacestr( spacestr_default ),
    _spacesep( spacesep_default ){}
};

/**
 * Calling the pure-virtual str() method of the format for output.
 */
inline std::ostream&
operator<<( std::ostream& os, const format& x )
{
  os << x.str();
  return os;
}

/**
 * @brief decimal representation of double.
 */
class decimal : public format
{
public:
  decimal(
    const double input,
    const int    p = -precision_default ) :
    _input( input )
  {
    precision( p );
  }

  virtual
  ~decimal (){}

  virtual std::string str() const;

private:
  double _input;
};

/**
 * @brief scientific notation of double.
 */
class scientific : public format
{
public:
  scientific( const double input, const unsigned p );
  virtual
  ~scientific  (){}

  virtual std::string str() const;

private:
  double _mant;
  int _exp;
};

/** @} */

}/* base */

}/* fmt */

}/* usr */

#endif/* end of include guard: USERUTILS_COMMON_STLUTILS_FORMAT_HPP */
