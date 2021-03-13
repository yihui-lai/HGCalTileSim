@addtogroup ostreamutils

This also includes classes for better outputting formats for the numbers
(double) that are more in-line with what is required for CMS publications.
All the could be altered is implemented in the usr::fmt::base::format class.
Including setting number of units to be displayed after the decimal point, and
additional spacing characters to be inserted every "n" digits.

The class is designed such that one can use it like:
```cpp
std::cout << decimal(mynumber).precision(4).spacesep(3).spacestr("..") ;
```
defaults for what the space separation string and number of digits could also
be set by the user via the global variables. For more details, see the
documentation in the implemented classes.


@class usr::fmt::base::format
@details

@class usr::fmt::base::decimal
@details
If the precision is set to a negative value, the absolute value of the given
precision will be used with the exception that the trailing zeros after the
decimal point will be stripped off.


@class usr::fmt::base::scientific
@details
Representing the number as the latex compatible string:

```
 "a \\times 10^{b}"
```

with \f$1 \leq |a| < 10\f$ and \f$b\neq0\f$. In the case that the calculated
exponent value is \f$0\f$, then output would be identical as that of using
the base::decimal class.
