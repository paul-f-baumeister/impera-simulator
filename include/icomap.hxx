#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <algorithm> // std::min

#include "data_view.hxx" // view2D<T>, view3D<T>
#include "inline_tools.hxx" // set, add_product
#include "status.hxx" // status_t

  typedef struct {
      uint32_t neighbors[6];
      float    resources;
      int8_t n_neighbors;
  } GridCell_t;

namespace icomap {

  inline constexpr
  size_t rhomb_edge(unsigned const Level) { return (1ul << Level); }
  
  inline constexpr
  size_t map_width(unsigned const Level) { return 10*rhomb_edge(Level); }
  
  inline constexpr 
  size_t map_height(unsigned const Level) { return 3*rhomb_edge(Level) + 1; }

  inline constexpr
  size_t ico_index(unsigned const Level, unsigned const i10, int const iSE, int const iNE) { 
      return iNE + rhomb_edge(Level)*(iSE + rhomb_edge(Level)*i10); }

  inline constexpr
  size_t n_ico_vertices(unsigned const Level) { return ico_index(Level, 10, 0, 2); } // 2 because we add NorthPole and SouthPole

  inline constexpr
  size_t pole_ico_index(unsigned const Level, char const sn) { // SouthPole: pole_ico_index(L,'S');
      return ico_index(Level, 10, 0, ('n' == (sn | 32))); }      // NorthPole: pole_ico_index(L,'N');


  template <typename real_t, typename real_in_t>
  status_t create_world_map(
        view3D<real_t> & map // data layout [Height][Width][stride]
      , view2D<real_in_t> const & dat // data layout [ico_index][stride]
      , int const Level
      , real_t const zero=0 // neutral element
      , real_t const f=1 // scale factor
      , bool const icosahedral=false // true: show map with corners from the ico-grid in it
  ) {
      // this map mimiques a Mollweide map
      size_t const Width = map_width(Level);
      size_t const tL = rhomb_edge(Level);

      size_t const stride = std::min(dat.stride(), map.stride());

      set(map.data(), map_height(Level)*Width*map.stride(), zero);
      
      auto const denom = 1./(2*tL - 1.);
      
      for (int i10 = 0; i10 < 10; ++i10) { // rhombs
          int const i01 = i10 & 0x1; // == i10 % 2; // 0: northern rhomb, 1: southern rhomb
// #ifdef MOLLWEIDE
//        int const im5 = (i10 - i01 + 5) % 10 - 5;
// #endif
          for (int iSE = 0; iSE < tL; ++iSE) { // south east direction
              for (int iNE = 0; iNE < tL; ++iNE) { // north east direction
                  auto const idx = ico_index(Level, i10, iSE, iNE);
                  int const iS = (i01 + 1)*tL + iSE - iNE,
                            iE = (i10 + 4)*tL + iSE + iNE; // the +4 takes care of the Europe-centered map with Greenwich in the very center
                            // or (if the map is shifted by -4 degrees) Madrid is right in the middle
// #ifdef MOLLWEIDE
//                   int const ipolar = (iSE - iNE)*(1 - 2*i01);
//                   if (ipolar < 0) { iE += im5*ipolar; }
// #endif
                  if (icosahedral || (((0 == i01) && (iSE > iNE)) || ((1 == i01) && (iSE < iNE)))) {
                      // determine 2 pixel per vertex to mimique the hex pattern, need to elongate 
                      // the picture by sqrt(3) in the y-direction to get an undistorted hex-map
                      int const iE0 = (iE + 0) % Width,
                                iE1 = (iE + 1) % Width;
                      set(map(iS,iE0), stride, dat[idx], f);
                      set(map(iS,iE1), stride, dat[idx], f);
                  } // central latitudes
              } // iNE
          } // iSE
          
          if (!icosahedral) {
              // interpolate the upper and lower triangles

              for (int iL = 1; iL <= tL; ++iL) { // towards the north/south pole
                  int const iS = i01 ? (3*tL - iL) : iL;
                  float const fE = (iL - 1)*denom; // (iL - 1)/(2*tL - 1.);
                  for (int iE = 0; iE < 2*tL; ++iE) { // east direction
                      float const EE = fE*iE;

                      int iSE0, iNE0, iSE1, iNE1;
                      float w0, w1;
                      if (i01) {
                          // southern rhomb
                          // linear interpolate between the points
                          //   dat[i10][tL - iS    ][0]     --> maps to iE=0
                          //   dat[i10][tL - iS + 1][1]
                          //   ...
                          //   dat[i10][tL - 1][iN - 1]     --> maps to iE=2*tL-1
                          iNE0 = int(EE), // floor
                          iNE1 = iNE0 + 1;
                          w1 = EE - iNE0;
                          w0 = 1 - w1;
                          iSE0 = tL - iL + iNE0,
                          iSE1 = tL - iL + iNE1;
                      } else {
                          // northern rhomb
                          // linear interpolate between the points
                          //   dat[i10][0][tL - iS    ]     --> maps to iE=0
                          //   dat[i10][1][tL - iS + 1]
                          //   ...
                          //   dat[i10][iS - 1][tL - 1]     --> maps to iE=2*tL-1
                          iSE0 = int(EE), // floor
                          iSE1 = iSE0 + 1;
                          w1 = EE - iSE0;
                          w0 = 1 - w1;
                          iNE0 = tL - iS + iSE0,
                          iNE1 = tL - iS + iSE1;
                      } // i01

                      // from here the same as in the block above
                      auto const idx0 = ico_index(Level, i10, iSE0, iNE0),
                                 idx1 = ico_index(Level, i10, iSE1, iNE1);
                      int const iEast = (iE + (i10 + 4)*tL) % Width;
                      set(map(iS,iEast), stride, dat[idx0], f*w0);
                      if (w1 > 0) add_product(map(iS,iEast), stride, dat[idx1], f*w1);
                  } // iE
              } // iL
          } // not icosahedral

      } // i10

      return 0;
  } // create_world_map
  
  
  status_t inline test_create_world_map(int const echo=1, int const Level=3) {
      status_t stat(0);
      view2D<char> input(n_ico_vertices(Level), 1, '@');
      for (size_t i = 0; i < n_ico_vertices(Level); ++i) {
          input(i,0) = 'a' + (i & 15);
      } // i
      view3D<float> output(map_height(Level), map_width(Level), 1);
      stat += create_world_map(output, input, Level, float(' '), 1.f, false);
      if (echo < 1) return stat; // early return
      std::printf("\n# Ico map %d x %d for Level %i\n", map_height(Level), map_width(Level), Level);
      for (int south = 0; south < map_height(Level); ++south) {
          std::printf("# ");
          for (int east = 0; east < map_width(Level); ++east) {
              std::printf("%c", char(output(south, east, 0)));
          } // east
          std::printf(" #\n");
      } // south

      std::printf("\n# Rectangular map %d x %d for Level %i\n", map_height(Level), map_width(Level), Level);
      view3D<float> square(map_height(Level), map_width(Level), 1);
      stat += create_world_map(square, input, Level, 0.f);
      for (int south = 0; south < map_height(Level); ++south) {
          std::printf("# ");
          for (int east = 0; east < map_width(Level); ++east) {
              std::printf("%c", char(square(south, east, 0)));
          } // east
          std::printf(" #\n");
      } // south

      return stat;
  } // test_create_world_map
  
  status_t inline all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_create_world_map(echo);
      return stat;
  } // all_tests

} // namespace icomap
