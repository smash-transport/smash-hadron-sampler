#include <array>

#ifndef INCLUDE_CONST_H_
#define INCLUDE_CONST_H_

#include <cmath>

// some general constants etc

constexpr double gevtofm = 5.067728853;
constexpr double hbarC = 1. / gevtofm;
constexpr double C_Feq = (std::pow(0.5 / M_PI / hbarC, 3));
constexpr double small_value = 1.e-10;
constexpr double gevtofm3_2pi2 = std::pow(gevtofm, 3) / (2. * M_PI * M_PI);

// Non zero components of the metric tensor
const std::array<int, 4> metric = {1, -1, -1, -1};

#endif // INCLUDE_CONST_H_
