#pragma once

#include "status.hxx" // status_t

namespace impera {

  status_t read_resource_file(std::vector<GridCell_t> & resources, char const *filename=nullptr, int const echo=0);

  inline std::vector<GridCell_t> read_resource_file(char const *filename=nullptr, int const echo=0) {
      std::vector<GridCell_t> resources;
      read_resource_file(resources, filename, echo);
      return resources;
  } // read_resource_file

  // interface for running the simulation in the background of a window GUI
  void* run_some_days(int const Nspecies, int const ndays=1, int const echo=9, int const fp_bits=64);
  int constexpr MemoryCleanup = -1; // passed into ndays

  // start the simulation in the terminal
  status_t all_tests(int const echo=1);

} // namespace impera
