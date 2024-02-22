#include "gen.h"

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
