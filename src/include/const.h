#include <array>

#ifndef INCLUDE_CONST_H_
#define INCLUDE_CONST_H_

// some general constants etc

const double gevtofm = 5.067728853;
const double hbarC = 1. / 5.067728853;
const double small_value = 1.e-10;

// 4d Levi-Civita tensor with indices in the range [0,3]
int epsilon_tensor(int i, int j, int k, int l);

// Non zero components of the metric tensor
const std::array<int, 4> metric = {1, -1, -1, -1};

#endif  // INCLUDE_CONST_H_
