#pragma once

namespace impera {

  status_t read_resource_file(std::vector<GridCell_t> & resources, char const *filename, int const echo=0);

  void* run_one_day(int const Nspecies, int const echo=9);

  int all_tests(int const echo=1);

} // namespace impera
