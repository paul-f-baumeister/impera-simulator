#pragma once

#include "status.hxx" // status_t

namespace window {

  status_t init(int argc, char *argv[]);

  void display3D(int Level=-1, float const *rgb=nullptr); // fill arrays for display

  status_t finalize();

  status_t all_tests(int const echo=0);

} // namespace window
