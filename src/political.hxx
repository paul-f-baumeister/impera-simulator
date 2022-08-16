#pragma once

#include <cmath> // std::sin, std::cos, std::abs
#include <algorithm> // std::min, std::max

#include "inline_tools.hxx" // pow2, set
#include "data_view.hxx" // view2D<T>
#include "color.hxx" // ::def, ::colorchar

namespace political {

  inline double check_political_colors(float const pcs[], unsigned const Nspecies) {
      double min_dist2{9e99};
      for (int is = 0; is < Nspecies; ++is) {
          for (int js = 0; js < is; ++js) { // compare to all previous
              double dist2{0};
              for (int rgb = 0; rgb < 3; ++rgb) {
                  dist2 += pow2(pcs[is*4 + rgb] - pcs[js*4 + rgb]);
              } // rgb
              std::min(min_dist2, dist2); // make sure the start colors are mutually different
          } // js
      } // is
      return min_dist2;
  } // check_political_colors
  
  inline double init_political_colors(float pcs[], unsigned const Nspecies) {
      for (int is = 0; is < Nspecies; ++is) {
          for (int rgb = 0; rgb < 3; ++rgb) {
              pcs[is*4 + rgb] = .875f*(rgb == (is % 3)) + .124f*std::sin(is) + .01f*std::cos(is); // init
          } // rgb
      } // is
      return check_political_colors(pcs, Nspecies);
  } // init_political_colors

  inline double update_political_colors(float pcs[], unsigned const Nspecies
            , double const overlap[] // assumed matrix shape [Nspecies][Nspecies], only off-diagonal elements are relevant
            , unsigned const Nsteps=99) {
      view2D<float> force(Nspecies, 4);
      float constexpr dt = 1e-3;
      float constexpr bw = 1e2; // control here how strong the force is that keeps the color saturation high
      float largest_displacement{1}, largest_force;
      unsigned step;
      for (step = 0; (step < Nsteps) && (largest_displacement > 2e-6); ++step) {
        
          // force maximing the color saturation
          for (int is = 0; is < Nspecies; ++is) {
              float const posi[] = {pcs[is*4], pcs[is*4 + 1], pcs[is*4 + 2]};
              auto const bw_force = bw*(1 - (posi[0] + posi[1] + posi[2])); // 1:slightly darker than grey
              set(force[is], 4, bw_force); // init force
          } // is

          // repulsive pair forces
          for (int is = 0; is < Nspecies; ++is) {
              float const posi[] = {pcs[is*4], pcs[is*4 + 1], pcs[is*4 + 2]};
              for (int js = 0; js < is; ++js) { // compare to all previous
                  float const posj[] = {pcs[js*4], pcs[js*4 + 1], pcs[js*4 + 2]};
                  float const diff[] = {posj[0] - posi[0], posj[1] - posi[1], posj[2] - posi[2]};
                  float const dist2 = pow2(diff[0]) + pow2(diff[1]) + pow2(diff[2]) + 0.01; // 0.01: some safety
                  float const ovl = (Nspecies > 6) ? overlap[is*Nspecies + js] : 0;
                  float const extra_force = (ovl == ovl) ? ovl : 0;
                  float const factor = 1.f/pow2(dist2) + extra_force;
                  add_product(force[is], 3, diff, -factor);
                  add_product(force[js], 3, diff,  factor);
              } // js
          } // is

          largest_force = 0;
          largest_displacement = 0;
          
          // relax
          for (int is = 0; is < Nspecies; ++is) {
              for (int rgb = 0; rgb < 3; ++rgb) {
                  float const oldpos = pcs[is*4 + rgb];
                  float const newpos = oldpos + dt*force(is,rgb); // suggested new position
                  pcs[is*4 + rgb] = std::min(std::max(0.f, newpos), 1.f);
                  largest_force = std::max(largest_force, std::abs(force(is,rgb)));
                  largest_displacement = std::max(largest_displacement, std::abs(pcs[is*4 + rgb] - oldpos));
              } // rgb
          } // is
          
          // this algorithm for Nspecies=3 converges to pure RGB, for Nspecies=6, magenta,cyan and yellow are added
          // so the color find themselfs in the 6 vertices of the RGB cube with are neither black (000) nor white (111)
          // for more than 6 species, we add the overlap to the repulsive force
          // printf("# in step #%i, largest force is %g and largest displacement is %g\n", step, largest_force, largest_displacement);
      } // step
      printf("# after %d steps, largest force is %g and largest displacement is %g\n", step, largest_force, largest_displacement);

      for (int is = 0; is < Nspecies; ++is) {
          auto const cc = color::colorchar(&pcs[is*4]);
          printf("# relaxed political color for species%3i %s       %s\n", is, &cc, color::def);
      } // is
      
      return check_political_colors(pcs, Nspecies);
  } // update_political_colors
  
  inline int all_tests(int const echo=1) { return -1; }

} // namespace political