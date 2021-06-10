@defgroup Common Common
@brief A collection of functions for simpler experience using STL and boost functions.

@details The [STL](http://en.cppreference.com/w/) and
[boost](https://www.boost.org/) template libraries are widely-available,
incredibly versatile library that could be used as a base for many programming
tasks. The two libraries being designed to be general purpose however makes using
the libraries extremely verbose when writing the code, and could be much
simplified if one already knows that your code would be using some specific use
cases.

This sub-packages aims to do just that, specializing the STL and boost libraries
to cases most commonly used during my HEP analysis trip, as well as providing
additional classes that are of more general purpose with only STL interactions.

Aside from STL and boost function enhancements, this library also provides some
ROOT container manipulator classes and functions that are handy for the latter
classes in the more specialized libraries.
