#pragma once

#include "status.hxx" // status_t

namespace window {
  
  status_t init(int argc=0, char *argv[]=nullptr);
  status_t finalize();
  
  void display3D(int Level=-1, float const *rgb=nullptr);

  status_t all_tests(int const echo=0);

} // namespace window
