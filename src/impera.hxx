#pragma once

namespace impera {

  status_t read_resource_file(std::vector<GridCell_t> & resources, char const *filename=nullptr, int const echo=0);

  inline std::vector<GridCell_t> read_resource_file(char const *filename=nullptr, int const echo=0) {
      std::vector<GridCell_t> resources;
      read_resource_file(resources, filename, echo);
      return resources;
  } // read_resource_file

  // interface for running the simulation in the background of a window GUI
  void* run_one_day(int const Nspecies, int const echo=9);

  // start the simulation in the terminal
  int all_tests(int const echo=1);

} // namespace impera
