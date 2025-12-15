#include <array>

#ifndef INCLUDE_CONST_H_
#define INCLUDE_CONST_H_

#include <cmath>

// some general constants etc

const double gevtofm = 5.067728853;
const double hbarC = 1. / gevtofm;
const double C_Feq = (pow(0.5 / M_PI / hbarC, 3));
const double small_value = 1.e-10;

// Non zero components of the metric tensor
const std::array<int, 4> metric = {1, -1, -1, -1};

#endif // INCLUDE_CONST_H_
