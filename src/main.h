// #error Please comment out the next two lines under linux, then comment this error
// #include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
#ifdef _WIN32
#include <windows.h> // Include if under windows
#endif

#include "quaternion.h"

/* prototype of the function you need to optimize */
template <typename T_return, class... T_args>
using compute_function = T_return (*)(T_args... args);
