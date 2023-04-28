#pragma once

#include <string>
#include "quaternion.h"

typedef void(*comp_func)(quaternion_t x[N],  quaternion_t y[N], quaternion_t A[N][N]);
void add_function(comp_func f, std::string name, int flop);
