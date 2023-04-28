// #error Please comment out the next two lines under linux, then comment this error
// #include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
// #include <windows.h> // Include if under windows

#include "quaternion.h"

/* prototype of the function you need to optimize */
typedef void(*comp_func)(quaternion_t x[N], quaternion_t y[N], quaternion_t A[N][N]);
