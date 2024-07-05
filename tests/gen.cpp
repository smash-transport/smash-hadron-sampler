#include "gen.h"

#include <ctime>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "virtest/vir/test.h"

using namespace gen;

TEST(index44) {
  // Call with all valid values and expect passing
  int indices[4][4] = {{0, 1, 3, 6}, {1, 2, 4, 7}, {3, 4, 5, 8}, {6, 7, 8, 9}};
  for (int row = 0; row < 4; row++) {
    for (int column = 0; column < 4; column++) {
      VERIFY(index44(row, column) == indices[row][column]);
    }
  }
}

TEST(index44_index_out_of_range) {
  int valid_index[4] = {0, 1, 2, 3};
  int invalid_index[4] = {-2, -1, 4, 5};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      // first index out of range
      try {
        index44(invalid_index[i], valid_index[j]);
        std::cout << "index44 unexpectedly passed with first argument "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
      // second index out of range
      try {
        index44(valid_index[i], invalid_index[j]);
        std::cout << "index44 unexpectedly passed with second argument "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
      // both indices out of range
      try {
        index44(invalid_index[i], invalid_index[j]);
        std::cout << "index44 unexpectedly passed with both arguments "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
    }
  }
}
