#ifndef VORTICITY_H
#define VORTICITY_H

#include <array>
#include <cstddef>
#include <memory>
#include <optional>

namespace gen {
struct element;
}

class Vorticity {
 public:
  // Default constructor nulls all vorticity components
  Vorticity() : vorticity_({0.0}) {}

  using reference = double&;
  using const_reference = const double&;

  // access the component i of the vorticity tensor as linear array
  reference operator[](size_t index) { return vorticity_[index]; }

  // Check that the vorticity file exists
  static void ensure_vorticity_file_exists_and_check_format();

  // Check that the freezeout is given in the extended format which contains
  // the energy density
  static void ensure_extended_freezeout_is_given();

  // Set the number of corona cells in the freezeout surface
  static void set_number_of_corona_cells();

  // const overload of the [] operator
  const_reference operator[](size_t index) const { return vorticity_[index]; }

  // Set vorticity_ to a new array
  void set_vorticity(const std::array<double, 16>& new_vorticity) {
    vorticity_ = new_vorticity;
  }

  // Get the vorticity tensor as a 1D array
  std::array<double, 16> get_vorticity() const { return vorticity_; }

  // given the complete freezeout surface, this function sets the vorticity
  // tensor in all surface cells from the vorticity file
  static void set_vorticity_in_surface_cells(gen::element* surf, int N);

  // return a component of the vorticity tensor as a 4x4 matrix
  reference at(int i, int j) { return vorticity_[i * 4 + j]; }
  const_reference at(int i, int j) const { return vorticity_[i * 4 + j]; }

  // Boost the vorticity tensor to the fluid rest frame with given boost matrix
  void boost_vorticity_to_fluid_rest_frame(const double (&boostMatrix)[4][4]);

 private:
  // number of corona cells in the freezeout surface. As the freezeout surface
  // starts with the corona cells, we need the number to set the vorticity
  // tensor in all corona cells to zero.
  static int num_corona_cells_;

  // components of the vorticity tensor
  std::array<double, 16> vorticity_;
};

#endif  // VORTICITY_H
