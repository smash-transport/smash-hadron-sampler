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

TEST(favored_spin_projection_valid) {
  MinMax limits = {-1.0, 1.0};
  int spin;

  // Spin 1/2 -> s = 1
  for (double vorticity = limits.minimum; vorticity <= limits.maximum;
       vorticity += 0.01) {
    spin = get_favored_spin_projection_in_cell(vorticity, limits, 1);
    if (vorticity < 0.0) {
      VERIFY(spin == -1);
    } else {
      VERIFY(spin == 1);
    }
  }
  // Value exactly at the binning edges
  spin = get_favored_spin_projection_in_cell(-1.0, limits, 1);
  VERIFY(spin == -1);
  spin = get_favored_spin_projection_in_cell(0.0, limits, 1);
  VERIFY(spin == 1);
  spin = get_favored_spin_projection_in_cell(1.0, limits, 1);
  VERIFY(spin == 1);

  // Spin 1 -> s = 2
  const double bin_width = (limits.maximum - limits.minimum) / 3.;

  for (double vorticity = limits.minimum; vorticity <= limits.maximum;
       vorticity += 0.01) {
    spin = get_favored_spin_projection_in_cell(vorticity, limits, 2);
    if (vorticity < (limits.minimum + bin_width)) {
      VERIFY(spin == -2);
    } else if ((limits.minimum + bin_width) <= vorticity &&
               vorticity < (limits.minimum + (2. * bin_width))) {
      VERIFY(spin == 0);
    } else if ((limits.minimum + (2. * bin_width)) <= vorticity &&
               vorticity < limits.maximum) {
      VERIFY(spin == 2);
    }
  }
}