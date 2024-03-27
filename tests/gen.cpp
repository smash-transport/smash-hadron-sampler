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

TEST(vorticity_validity) {
  // Unset vorticity extrema
  MinMax limits_1 = {-5.9, std::numeric_limits<double>::quiet_NaN()};
  MinMax limits_2 = {std::numeric_limits<double>::quiet_NaN(), 5.9};
  MinMax limits_3 = {std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN()};
  MinMax limits_array[3] = {limits_1, limits_2, limits_3};

  for (MinMax limits : limits_array) {
    try {
      check_if_vorticity_values_are_valid(1.9, limits);
      std::cout << "Value Error: check_if_vorticity_values_are_valid() "
                << " unexpectedly passed with invalid vorticity extrema "
                << std::endl;
      FAIL();
    } catch (std::invalid_argument &e) {
      // Exception was caught as expected
    }
  }
  // Vorticity in cell exceeds vorticity extrema
  MinMax limits = {-1.0, 1.0};
  double invalid_vorticity[6] = {-999.9, -2.5, -1.01, 1.01, 2.5, 999.9};
  for (double vorticity : invalid_vorticity) {
    try {
      check_if_vorticity_values_are_valid(vorticity, limits);
      std::cout
          << "Initialization Error:  get_favored_spin_projection_in_cell() "
          << " unexpectedly passed with invalid vorticity extrema "
          << std::endl;
      FAIL();
    } catch (std::invalid_argument &e) {
      // Exception was caught as expected
    }
  }
}

TEST(favored_spin_projection) {
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

TEST(update_vorticity_extrema) {
  const int num_values = 1000;
  MinMax vorticity_extrema = {std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN()};
  std::vector<double> values(num_values);

  // Fill the vector with random values (negative to positive)
  std::default_random_engine generator(std::time(nullptr));
  std::uniform_real_distribution<double> distribution(-10.0, 10.0);
  for (int i = 0; i < num_values; ++i) {
    values[i] = distribution(generator);
    update_vorticity_extrema(i, values[i], vorticity_extrema);
  }
  // Find both minimum and maximum values of the vector
  auto min_max_pair = std::minmax_element(values.begin(), values.end());
  const double min_value = *min_max_pair.first;
  const double max_value = *min_max_pair.second;

  VERIFY(vorticity_extrema.minimum == min_value);
  VERIFY(vorticity_extrema.maximum == max_value);
}

TEST(spin_sampling_return_values) {
  const int num_samples = 50000;
  // Spin 0
  for (int i = 0; i < num_samples; i++) {
    int state = sample_spin_projection(0, 0, 0.0);
    VERIFY(state == 0);
  }
  // Spin 1
  for (int i = 0; i < num_samples; i++) {
    int state = sample_spin_projection(1, 1, 0.2);
    VERIFY(state == -1 || state == 1);
  }
  // Spin 2
  for (int i = 0; i < num_samples; i++) {
    int state = sample_spin_projection(2, -2, 0.2);
    VERIFY(state == -2 || state == 0 || state == 2);
  }
  // Spin 3
  for (int i = 0; i < num_samples; i++) {
    int state = sample_spin_projection(3, 3, 0.2);
    VERIFY(state == -3 || state == -1 || state == 1 || state == 3);
  }
  // Spin 4
  for (int i = 0; i < num_samples; i++) {
    int state = sample_spin_projection(4, -4, 0.2);
    VERIFY(state == -4 || state == -2 || state == 0 || state == 2 ||
           state == 4);
  }
  // Spin 5
  for (int i = 0; i < num_samples; i++) {
    int state = sample_spin_projection(5, 3, 0.2);
    VERIFY(state == -5 || state == -3 || state == -1 || state == 1 ||
           state == 3 || state == 5);
  }
}

// Given a vector of num_samples counts distributed among its elements, check
// that all counts but those in the element `skip_index`coincide with
// `probability` up to a given `precision` in percent.
static bool check_other_elements_than_given_index(
    const std::vector<int> &values, int skip_index, double probability,
    int num_samples, double precision) {
  // No other element to check
  if (values.size() == 1) {
    return true;
  } else {
    size_t index_to_skip = static_cast<size_t>(skip_index);

    // Check if index_to_skip is within the bounds of the vector
    if (index_to_skip >= values.size() || index_to_skip < 0) {
      throw std::out_of_range("Index to skip is out of range.");
    }

    // Iterate over all elements except values[favored_spin]
    for (size_t i = 0; i < values.size(); ++i) {
      if (i != index_to_skip) {
        double disfavored_probability =
            static_cast<double>(values[i]) / num_samples;
        if (std::abs(disfavored_probability - probability) > precision) {
          return false;
        }
      }
    }
    return true;
  }
}

TEST(spin_sampling_probability) {
  const double precision = 0.01;
  const double epsilon = 0.0001;
  const int num_samples = 600000;
  const double polarization = 0.2;

  // Loop through different spin states
  for (int spin = 0; spin < 6; spin++) {
    double num_states = spin + 1;
    double base_probability = 1.0 / num_states;

    // Calculate the probabilities for the favored state and all others
    double favored_probability;
    double disfavored_probability;
    if (spin == 0) {
      favored_probability = 1.0;
      disfavored_probability = 0.0;
    } else {
      favored_probability = (1. + polarization) * base_probability;
      disfavored_probability = (1.0 - favored_probability) / (num_states - 1);
    }

    // Create and fill the vector with allowed states (even values from -spin to
    // spin)
    std::vector<int> valid_states(num_states);
    for (int i = 0; i < num_states; ++i) {
      valid_states[i] = -spin + 2 * i;
    }

    // Create a vector to store the counts for each allowed state
    std::vector<int> counts(spin + 1, 0);

    for (int favored_spin = 0; favored_spin <= spin; favored_spin++) {
      // Loop over samples to generate statistics and fill counts
      for (int sample = 0; sample < num_samples; sample++) {
        int state = sample_spin_projection(spin, valid_states[favored_spin],
                                           polarization);
        // Take into account step-size of 2 between states
        int index = (state + spin) / 2;
        counts[index]++;
      }

      // Check that sum of all counts coincides with num_samples
      int sum_counts = std::accumulate(counts.begin(), counts.end(), 0);
      VERIFY(sum_counts == num_samples);

      // Check that the favored probability is sampled correctly
      double probability =
          static_cast<double>(counts[favored_spin]) / num_samples;
      VERIFY(std::abs(probability - favored_probability) < precision);

      // Check that all other probabilities coincide with the disfavored
      // probability
      bool are_other_probabilities_disfavored =
          check_other_elements_than_given_index(counts, favored_spin,
                                                disfavored_probability,
                                                num_samples, precision);
      VERIFY(are_other_probabilities_disfavored);

      // Check that the sum of all probabilities is normalized to 1
      double total_probability = 0.0;
      for (int count : counts) {
        total_probability += static_cast<double>(count) / num_samples;
      }
      VERIFY(std::abs(total_probability - 1.0) <= epsilon);

      // Reset counts for next favored spin index
      counts.assign(counts.size(), 0);
    }
  }
}

TEST(spin_sampling_invalid_initialization) {
  // Invalid spin
  for (int i = -1; i > -20; i--) {
    try {
      sample_spin_projection(-1, -1, 0.1);
      std::cout << "Value Error: sample_spin_projection() unexpectedly passed "
                << "with negative spin" << std::endl;
    } catch (std::invalid_argument &e) {
      // Exception was caught as expected
    }
  }

  // favored_spin_projection out of range
  for (int spin = 0; spin < 10; spin++) {
    for (int favored_spin = spin + 1; favored_spin < 20; favored_spin++) {
      try {
        sample_spin_projection(spin, favored_spin, 0.1);
        std::cout
            << "Value Error: sample_spin_projection() unexpectedly passed "
            << "with favored_spin_projection out of range" << std::endl;
      } catch (std::invalid_argument &e) {
        // Exception was caught as expected
      }
    }
  }

  // polarization_percentage larger than overall probability for spin > 0
  for (int spin = 1; spin < 7; spin++) {
    double max_polarization = static_cast<double>(spin);
    for (double polarization = max_polarization + 0.1;
         polarization <= max_polarization + 1.0; polarization += 0.1) {
      try {
        sample_spin_projection(spin, spin, polarization);
        std::cout
            << "Value Error: sample_spin_projection() unexpectedly passed "
            << "with polarization exceeding the maximum possible" << std::endl;
      } catch (std::invalid_argument &e) {
        // Exception was caught as expected
      }
    }
  }
}
