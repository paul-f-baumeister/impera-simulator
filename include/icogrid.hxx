#pragma once

#include <cstdio> // std::printf, ::fopen, ::fprintf, ::fclose
#include <cassert> // assert
#include <cmath> // std::cos, ::sin, M_PI, ::sqrt
#include <ctime> // std::time, ::time_t, ::put_time
#include <iomanip> // std::put_time
#include <vector> // std::vector<T>

#include "data_view.hxx" // view2D<T>, view3D<T>
#include "inline_tools.hxx" // pow2
#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "warnings.hxx" // warn

namespace icogrid {
  // Icosahedral grid construction

  double constexpr Degrees72 = 2*M_PI/5.; // angle between two pentagon points

  inline void base_point(double v[3], int const i10) {
      double const factor = std::sqrt(4/5.);
      double const angle = (i10 - 1)*(M_PI/5.); // the -1 makes the 0-th meridian pass in the middle through the 0-th rhomb (better to represent Europe)
//    double const angle = (i10 - 1 - 4./36.)*(M_PI/5.); // in order to fit SouthAmerica intirely into rhomb#3, shift by -4 degrees
      v[0] = factor*std::cos(angle);
      v[1] = factor*std::sin(angle);
      v[2] = factor*(0.5 - (i10 & 0x1)); // i & 0x1 is eqivalent to (i + 2)%2
      /* BasePoints
       *
       *    N   N   N   N   NorthPole
       *   / \ / \ / \ / \ / \
       *  0 - 2 - 4 - 6 - 8 - 0
       *   \ / \ / \ / \ / \ / \
       *    1 - 3 - 5 - 7 - 9 - 1
       *     \ / \ / \ / \ / \ /
       *      S   S   S   S   SouthPole
       */
  } // base_point of a rhomb in a doedecahedron inside a unit sphere

  inline constexpr
  size_t rhomb_edge(unsigned const Level) { return (1ul << Level); }

  inline constexpr
  size_t ico_index(unsigned const Level, unsigned const i10, unsigned const iSE, unsigned const iNE) {
      return (i10*rhomb_edge(Level) + iSE)*rhomb_edge(Level) + iNE; }


  inline constexpr
  size_t n_ico_triangles(unsigned const Level) { return 20*pow2(rhomb_edge(Level)); }

  inline constexpr
  size_t n_ico_vertices(unsigned const Level) { return ico_index(Level, 10, 0, 2); } // 2 because we add NorthPole and SouthPole

  inline constexpr
  size_t pole_ico_index(unsigned const Level, char const sn) { // NorthPole: pole_ico_index(L,'N');
      return  ico_index(Level, 10, 0, ('s' == (sn | 32))); }   // SouthPole: pole_ico_index(L,'S');


  inline size_t ico_index_wrap(unsigned const Level, unsigned const i10, unsigned const iSE, unsigned const iNE) {
      // Given we add 1 to the iSE or the iNE index of a grid point. Then we need to invoke this wrap function to
      // correct for the boundaries in case the point has left the rhomb.
 
      // std::printf("# ico_index_wrap(Level=%d, i10=%d, iSE=%d, iNE=%d)\n", Level, i10, iSE, iNE);
      int const tL = rhomb_edge(Level);
      /* RhombIndex (i10)
       *
       *    N   N   N   N   NorthPole
       *   / \ / \ / \ / \ / \
       *  x 0 x 2 x 4 x 6 x 8 x
       *   \ / \ / \ / \ / \ / \
       *    x 1 x 3 x 5 x 7 x 9 x
       *     \ / \ / \ / \ / \ /
       *      S   S   S   S   SouthPole
       */
      int j10, jSE, jNE;
      if (i10 & 0x1) {
          // i10 is odd, southern rhomb
          if (iSE == tL) {
              if (0 == iNE) return pole_ico_index(Level, 'S'); // south pole, i10=10, iSE=0, iNE=1
              // leave a southern rhomb to the east
              j10 = (i10 + 2) % 10; // hit another southern rhomb
              jNE = iSE - tL; // which is 0
              jSE = tL - iNE;
          } else if (iNE == tL) {
              // leave a southern rhomb to the north east
              j10 = (i10 + 1) % 10; // hit a northern rhomb
              jNE = iNE - tL;
              jSE = iSE;
          }
      } else {
          // i10 is even, northern rhomb
          if (iNE == tL) {
              if (0 == iSE) return pole_ico_index(Level, 'N'); // north pole, i10=10, iSE=0, iNE=0
              // leave a northern rhomb to the east
              j10 = (i10 + 2) % 10; // hit another northern rhomb
              jNE = tL - iSE;
              jSE = tL - iNE;
          } else if (iSE == tL) {
              // leave a northern rhomb the the south east
              j10 = (i10 + 1) % 10; // hit a southern rhomb
              jSE = iSE - tL;
              jNE = iNE;
          }
      }
      // std::printf("# ico_index_wrap(Level=%d, i10=%d, iSE=%d, iNE=%d) --> (%d, %d, %d)\n", Level, i10, iSE, iNE, j10, jSE, jNE);
      assert(0 <= j10); assert(j10 < 10);
      assert(0 <= jSE); assert(jSE < tL);
      assert(0 <= jNE); assert(jNE < tL);
      return (j10*tL + jSE)*tL + jNE; // == ico_index(j10, jSE, jNE);
  } // ico_index_wrap

  inline double norm2(double const x, double const y, double const z) { return pow2(x) + pow2(y) + pow2(z); }
  inline double norm2(double const v[3]) { return norm2(v[0], v[1], v[2]); }
  inline double normalize(double const v[3]) { return 1./std::sqrt(norm2(v)); }

  template <typename real_t> inline
  void crossProduct(real_t axb[3], real_t const a[3], real_t const b[3]) { 
      axb[0] = a[1] * b[2] - a[2] * b[1];
      axb[1] = a[2] * b[0] - a[0] * b[2];
      axb[2] = a[0] * b[1] - a[1] * b[0];
  } // crossProduct

//   inline int log2(size_t x) { int l2{-1}; while (x) { ++l2; x >>= 1; } return l2; } 

  template <typename real_t> inline
  status_t find_rhomb(double const lat, double const lon, view4D<real_t> const & vtx, int const echo=9) {
      double const degrees = 180/M_PI;
      double const critical_latitude = std::asin(std::sqrt(1./5.));
      double const deg36 = M_PI/5.; // 36 degrees
      printf("# latitude = %.1f longitude = %.1f degrees\n", lat*degrees, lon*degrees);
      double const l10 = lon/deg36;

      double const pos[3] = {std::cos(lon)*std::cos(lat), std::sin(lon)*std::cos(lat), std::sin(lat)};
      int j10{-1};
      if (std::abs(lat) > critical_latitude) {
          int const ins = (lat < 0); // 0: north, 1: south
          int const i10 = int(std::floor(l10) - ins + 10) % 10; // may need some shift
          j10 = i10;
      } else {
          int const i10 = int(std::floor(l10));
          // contruct a plane through points #i10 #(i10 + 1)%10 and the earth center
          // if the distance is positive, its rhomb #(i10 - 1)%10, else rhomb #i10
          double plane[3], v0[3], v1[3]; // normal vector
          base_point(v0, i10 + 0);
          base_point(v1, i10 + 1);
          crossProduct(plane, v0, v1);
          scale(plane, 3, normalize(plane));
          double const distance = dot_product(3, plane, pos);
          int const ins = i10 % 2;
          j10 = i10 + ((distance < 0) ^ ins);
      }

      int const tL = vtx.dim1() - 1;
      // now find the position inside the rhomb j10
      int iSE_start{0}, iSE_end{tL}, iSE = tL >> 1;
      int iNE_start{0}, iNE_end{tL}, iNE{-1};
      for (int iter = 0; iter < 4; ++iter) {

          // minimize in iNE
          {
              double d2min{9e9};
              for(int jNE = iNE_start; jNE <= iNE_end; ++jNE) {
                  double const *const v = vtx(j10, iSE, jNE);
                  double const dist2 = norm2(pos[0] - v[0], pos[1] - v[1], pos[2] - v[2]);
                  if (dist2 < d2min) { iNE = jNE; d2min = dist2; }
              } // jNE
              if (echo > 7) printf("# %i-th minimize found iNE=%i, distance^2= %g\n", iter, iNE, d2min);
          }

          // minimize in iSE
          {
              double d2min{9e9};
              for(int jSE = iSE_start; jSE <= iSE_end; ++jSE) {
                  double const *const v = vtx(j10, jSE, iNE);
                  double const dist2 = norm2(pos[0] - v[0], pos[1] - v[1], pos[2] - v[2]);
                  if (dist2 < d2min) { iSE = jSE; d2min = dist2; }
              } // jSE
              if (echo > 7) printf("# %i-th minimize found iSE=%i, distance^2= %g\n", iter, iSE, d2min);
          }

      } // iter

      return 0;
  } // find_rhomb



  inline int32_t global_coordinates(unsigned Level, unsigned i10, unsigned iSE, unsigned iNE) {
      assert(Level < 15); // 10*4^15 > 2^32
      // Level=14 --> 5,368,709,120 triangles, resolution ~308m
      // 2^32 ~= 4.3G --> uint32_t can index all 2,684,354,560 vertices
      //                  using int32_t, some indices will be negative!

//       int32_t constexpr id_north_pole = (0xa << 28);
//       int32_t constexpr id_south_pole = (0xa << 28) | 0x1;
      int32_t id = (i10 % 10) << 28;
      unsigned const edge = rhomb_edge(Level);
      assert(iSE < edge);
      assert(iNE < edge);
      // now assume that iSE and iNE are inside [0, 2^Level)
//       bool const southern_hemisphere = i10 & 0x1; // true for odd 
//       if (iSE == edge && southern_hemisphere)
      for(int L = 0; L < Level; ++L) {
          id |= (((iNE >> L) & 0x1) << (2*(14 - L)    ));
          id |= (((iSE >> L) & 0x1) << (2*(14 - L) + 1));
      } // L
      return id;
  } // global_coordinates

  template <typename int_t=int16_t>
  void global_coordinates(int_t i10_iSE_iNE[3], int32_t const id) {
      unsigned iNE{0}, iSE{0};
      for(int L = 0; L < 15; ++L) {
          iSE |= ((id >> (2*L + 1)) & 0x1) << L;
          iNE |= ((id >> (2*L    )) & 0x1) << L;
      } // L
      i10_iSE_iNE[0] = id >> 28; // i10
      i10_iSE_iNE[1] = iSE;
      i10_iSE_iNE[2] = iNE;
  } // global_coordinates

  template <typename real_t> inline
  status_t find_rhomb(double const xyz[3], view4D<real_t> const & vtx) {
      double const n2 = norm2(xyz);
      if (n2 < .5) return -1; // not even close to the unit sphere
      double const f = 1./std::sqrt(n2); // normalization factor
      double const v[3] = {xyz[0]*f, xyz[1]*f, xyz[2]*f};
      // now v is on the unit sphere
      double const lat = std::asin(v[2]), lon = std::atan2(v[1], v[0]);
      return find_rhomb(lat, lon, vtx);
  } // find_rhomb


  template <typename real_t>
  status_t inline generate(view4D<real_t> & vtx, unsigned const Level, int const echo=0) {
      assert(vtx.stride() >= 3);
      set(vtx, 10, real_t(0));

      if (0 == Level) {
          // for Level 0, we generate the 12 pentagon points of a icosahedron
          // which has 20 triangle faces

          double vtx12[12][3];
          for (int i10 = 0; i10 < 10; ++i10) {
              auto v = vtx12[i10];
              base_point(v, i10);
              double const norm2 = pow2(v[0]) + pow2(v[1]) + pow2(v[2]) - 1;
              if (echo > 17) printf("# %s i10=%i %.1e\n", __func__, i10, norm2);
              if (echo > 7) printf("# %s i10=%i %16.9f%16.9f%16.9f\n", __func__, i10, v[0], v[1], v[2]);
          } // i10

          for (char sn = 'N'; sn <= 'S'; sn += 'S' - 'N') {
              auto v = vtx12[pole_ico_index(0, sn)];
              v[0] = 0; v[1] = 0; v[2] = ('N' == sn)? 1 : -1;
          } // sn

          // expand 12 unique points into 10*2*2 points
          assert(2 == vtx.dim1());
          assert(2 == vtx.dim2());
          int const iNP = pole_ico_index(0, 'N');
          int const iSP = pole_ico_index(0, 'S');
          for (int i5 = 0; i5 < 5; ++i5) {
              {   // Northern rhomb (even)
                  unsigned const i10 = 2*i5 + 0;
                  set(vtx(i10, 0, 0), 3, vtx12[i10]);
                  set(vtx(i10, 0, 1), 3, vtx12[iNP]);
                  set(vtx(i10, 1, 0), 3, vtx12[(i10 + 1) % 10]);
                  set(vtx(i10, 1, 1), 3, vtx12[(i10 + 2) % 10]);
              }
              {   // Southern rhomb (odd)
                  unsigned const i10 = 2*i5 + 1;
                  set(vtx(i10, 0, 0), 3, vtx12[i10]);
                  set(vtx(i10, 0, 1), 3, vtx12[(i10 + 1) % 10]);
                  set(vtx(i10, 1, 0), 3, vtx12[iSP]);
                  set(vtx(i10, 1, 1), 3, vtx12[(i10 + 2) % 10]);
              }
          } // i5

          return 0;
      } // 0 == Level


      size_t const tL = rhomb_edge(Level - 1);
      view4D<real_t> vtx_prev(10, tL + 1, tL + 1, 3);
      // recursive invokation
      auto stat = generate(vtx_prev, Level - 1, echo - 1);

      assert(2*tL + 1 == vtx.dim1());
      assert(2*tL + 1 == vtx.dim2());

      double longest2{0}, shortest2{9e9}; // stats for the distances^2
      size_t nerrors{0};
      for (int i10 = 0; i10 < 10; ++i10) { // rhombs

          // copy the matching points (even iNE, even iSE) from the previous grid level
          for (int iSE = 0; iSE <= tL; ++iSE) { // south east direction
              for (int iNE = 0; iNE <= tL; ++iNE) { // north east direction
                  set(vtx(i10, 2*iSE, 2*iNE), 3, vtx_prev(i10, iSE, iNE));
              } // iNE
          } // iSE

          //=================================================================================
          // generate the new points in between:
          //=================================================================================

          real_t const one = 1;

          // create points (odd iNE, even iSE) from (iNE-1, iSE) and (iNE+1, iSE)
          for (int iSE = 0; iSE <= 2*tL; iSE += 2) { // south east direction
              for (int iNE = 1; iNE <= 2*tL; iNE += 2) { // north east direction
                  double v[3]; set(v, 3, vtx(i10, iSE, iNE - 1));
                  add_product(     v, 3, vtx(i10, iSE, iNE + 1), one);
                  set(vtx(i10, iSE, iNE), 3, v, normalize(v));
              } // iNE
          } // iSE

          // create points (even iNE, odd iSE) from (iNE, iSE-1) and (iNE, iSE+1)
          for (int iSE = 1; iSE <= 2*tL; iSE += 2) { // south east direction
              for (int iNE = 0; iNE <= 2*tL; iNE += 2) { // north east direction
                  double v[3]; set(v, 3, vtx(i10, iSE - 1, iNE));
                  add_product(     v, 3, vtx(i10, iSE + 1, iNE), one);
                  set(vtx(i10, iSE, iNE), 3, v, normalize(v));
              } // iNE
          } // iSE

          // create points (odd iNE, odd iSE) from (iNE-1, iSE-1) and (iNE+1, iSE+1)
          for (int iSE = 1; iSE <= 2*tL; iSE += 2) { // south east direction
              for (int iNE = 1; iNE <= 2*tL; iNE += 2) { // north east direction
                  double v[3]; set(v, 3, vtx(i10, iSE - 1, iNE - 1));
                  add_product(     v, 3, vtx(i10, iSE + 1, iNE + 1), one);
                  set(vtx(i10, iSE, iNE), 3, v, normalize(v));
              } // iNE
          } // iSE

          // DEBUG : count vectors which are not normalized
          for (int iSE = 0; iSE <= 2*tL; ++iSE) { // south east direction
              for (int iNE = 0; iNE <= 2*tL; ++iNE) { // north east direction
                  nerrors += (std::abs(norm2(vtx(i10, iSE, iNE)) - 1) > 1e-6);
              } // iNE
          } // iSE

          // distance stats
          for (int iSE = 0; iSE < 2*tL; ++iSE) { // south east direction
              for (int iNE = 0; iNE < 2*tL; ++iNE) { // north east direction
                  real_t dv[3];
                  {
                      set(dv, 3, vtx(i10, iSE, iNE));
                      add_product(dv, 3, vtx(i10, iSE + 0, iNE + 1), real_t(-1));
                      longest2  = std::max(longest2,  norm2(dv));
                      shortest2 = std::min(shortest2, norm2(dv));
                  }
                  {
                      set(dv, 3, vtx(i10, iSE, iNE));
                      add_product(dv, 3, vtx(i10, iSE + 1, iNE + 0), real_t(-1));
                      longest2  = std::max(longest2,  norm2(dv));
                      shortest2 = std::min(shortest2, norm2(dv));
                  }
                  {
                      set(dv, 3, vtx(i10, iSE, iNE));
                      add_product(dv, 3, vtx(i10, iSE + 1, iNE + 1), real_t(-1));
                      longest2  = std::max(longest2,  norm2(dv));
                      shortest2 = std::min(shortest2, norm2(dv));
                  }

              } // iNE
          } // iSE


      } // i10
      if (nerrors) printf("# %s Level=%i found %ld errors!\n", __func__, Level, nerrors);
      // else printf("# %s Level=%i found no errors!\n", __func__, Level); // DEBUG

      double const Radius = control::get("earth.radius", 6.3781e6); // radius of the Earth in meters
      auto const SegmentDistance = Radius/std::sin(Degrees72); // space distance between two pentagon centers
      if (echo > 10) {

          // show some characteristic number of this grid level
          size_t const nTriangles = 20 * (1ull << (2*Level)); // 20*4^L
          size_t const nEdges     = 12 * (1ull << (2*Level)); // 12*4^L // not sure
          auto const Circumference = 2*M_PI*Radius;
          auto const SurfaceArea = 2*Circumference*Radius;
          auto const Volume = SurfaceArea*Radius/3.;
          auto const BowDistance = Radius*2*std::asin(0.5/std::sin(Degrees72)); // surface distance on great circle between two pentagon centers

          if (1 == Level) {
              printf("# Earth Radius %.6e m, Circumference %.6e m, Surface Area %.6e m^2, Volume %.6e m^3\n",
                        Radius, Circumference, SurfaceArea, Volume);
              printf("# Segment Distance %.6e m, Bow Distance %.6e m\n", SegmentDistance, BowDistance);
          } // Level 1

          auto const TriangleArea = SurfaceArea/nTriangles;
          printf("# Grid-Level %2d:\n", Level);
          printf("#\tNumber of Triangles %.6f M, Edges %.6f M\n", nTriangles*1e-6, nEdges*1e-6);
          printf("#\tMemory Consumption %.6f GByte/float\n", nTriangles*4e-9);
          printf("#\tAverage Area of each Triangle %.3f m^2\n", TriangleArea);
          printf("#\tLength Scale of the Resolution %.3f m\n", std::sqrt(TriangleArea));
          printf("#\tApproximate Triangle Edge Length %.3f m\n", std::sqrt(4*TriangleArea/std::sqrt(3.)));
          printf("#\tApproximate Point Distance %.3f m\n", BowDistance/(1ull << Level));
          printf("#\tLongest and Shortest Distance %.3f and %.3f m\n", Radius*std::sqrt(longest2), Radius*std::sqrt(shortest2));
          printf("#\tLongest over Shortest ratio %.9f\n", std::sqrt(longest2 / shortest2));
          printf("#\n\n");
      } // echo


      auto const svg_filename = control::get("icogrid.export.svg.file", "");
      if (svg_filename && *svg_filename) {
          auto const svgf = std::fopen(svg_filename, "w");
          std::fprintf(svgf, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
          std::fprintf(svgf, "<!-- Creator: %s -->\n", __FILE__);
          auto const time_now = std::time(nullptr);
          std::fprintf(svgf, "<!-- CreationDate: %s -->\n", std::put_time(std::gmtime(&time_now), "%c %Z"));
// <!-- Magnification: 1 -->
          double const unit = SegmentDistance*1e-3; // in km, usually 6706.331 km
          std::fprintf(svgf, "<svg width=\"%dpt\" height=\"%dpt\" viewBox=\"%d %d %d %d\" "
              "xmlns=\"http://www.w3.org/2000/svg\">\n", 512, 256, 0, 0, int(6*unit), int(3*unit));

          double const East2_unit  =           0.5/tL*unit;
          double const North_unit = std::sqrt(.75)/tL*unit;
          // the eastward extend if a triangle is 2 East2_unit
          //   while the northward extend is only 1 North_unit

          for (int i10 = 0; i10 < 10; ++i10) {
              double const rhomb_East  = (0.25 + i10)     *tL*East2_unit;
              double const rhomb_North = (0.5 - (i10 % 2))*tL*North_unit;

              { // draw the outer triangles (fold and cut lines)
                  int8_t const from_to[5][2][2] = {{{0, 0}, {2,  0}},   // base point i10 --> base point i10 + 2
                                                   {{0, 0}, {1,  1}},   // base point i10 --> northern rhomb point
                                                   {{0, 0}, {1, -1}},   // base point i10 --> southern rhomb point
                                                   {{2, 0}, {1,  1}},   // base point i10 + 2 --> northern rhomb point
                                                   {{2, 0}, {1, -1}}};  // base point i10 + 2 --> southern rhomb point
                  for (int ft = 0; ft < 5; ++ft) {
                      auto const fr = from_to[ft][0];
                      auto const to = from_to[ft][1];
                      double const line[2][2] = {{fr[0]*tL*East2_unit + rhomb_East, fr[1]*tL*North_unit + rhomb_North},
                                                 {to[0]*tL*East2_unit + rhomb_East, to[1]*tL*North_unit + rhomb_North}};
                       std::fprintf(svgf, "<line stroke=\"black\" stroke-width=\"5\" stroke-linecap=\"round\" "
                                          "x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" />\n",
                                           line[0][0], line[0][1], line[1][0], line[1][1]);
                  } // ft
              } // black lines

              { // draw the inner triangles (mesh lines)

                  for (int iSE = 0; iSE < 2*tL; ++iSE) { // south east direction
                      for (int iNE = 0; iNE < 2*tL; ++iNE) { // north east direction
                          for (int ft = 0; ft < 2; ++ft) {

                          } // ft
                      } // iNE
                  } // iSE

              } // grey lines

          } // i10



          std::fprintf(svgf, "</svg>\n");
          stat += std::fclose(svgf);
      } // svg_filename

      return stat;
  } // generate

  status_t inline test_triangle(int const echo=1, unsigned const Level=3) {
      // find the symmetry multiplicity of equivalent points in a Level-subdivided triangle
      // exploiting 3-fold rotation about the center and 2-fold mirror planes

      int const edge = rhomb_edge(Level);

      view2D<int8_t> level(edge + 1, edge + 1, -1); // init levels as impossible
      for(int L = Level; L >= 0; --L) {
          int const inc = 1 << (Level - L); // incremenet
//        printf("# %s: Level %d of %d has increment %d\n", __func__, L, Level, inc);
          for(int i = 0; i <= edge; i += inc) {
              for(int j = 0; j <= i; j += inc) { // triangular loop
                  level(i,j) = L; // this grid point belongs to level L (and any higher level, as they are inclusive)
              } // j
          } // i
      } // L

      view2D<uint8_t> multi(edge + 1, edge + 1, 0); // init multiplicities
      // how to index points in that triangle? --> with a triangular loop
      for(int i = 0; i < edge; ++i) {
          for(int j = 0; j < i; ++j) { // triangular loop, lower triangular matrix
              multi(i,j) = 20; // there are 20 triangle faces in an icosahedron
          } // j
          multi(i,0) = 30; // edge owned by this triangle, there are 30 edges in an icosahedron
      } // i
      multi(0,0) = 12; // special point: the base point


      // rotating operation maps:
      //    00 --> ee
      //    e0 --> 00
      //    ee --> e0

      // rotation (0, 0) --> (4, 4)
      // rotation (1, 0) --> (3, 3)
      // rotation (1, 1) --> (4, 3)
      // rotation (2, 0) --> (2, 2)
      // rotation (2, 1) --> (3, 2)
      // rotation (2, 2) --> (4, 2)
      // rotation (3, 0) --> (1, 1)
      // rotation (3, 1) --> (2, 1)
      // rotation (3, 2) --> (3, 1)
      // rotation (3, 3) --> (4, 1)
      // rotation (4, 0) --> (0, 0)
      // rotation (4, 1) --> (1, 0)
      // rotation (4, 2) --> (2, 0)
      // rotation (4, 3) --> (3, 0)
      // rotation (4, 4) --> (4, 0)

      // implemented as (i,j) --> (e+j-i,e-i)

      // base vectors: (1,0) for i and (-.5, s34) for j
      // rotate by 120 degrees == 2*pi/3
      // using
      //    | -.5  s34 |
      //    | -s34 -.5 |
      // rotation about the center: 
      // the center is the center of mass between the 3 vertex points
      // 00, e0 and ee, i.e. 2e/3*(1,0) + e/3*(-.5, s34)
      // == e*(.5, s34/3)

#ifdef GEOMETRIC
      double const s34 = std::sqrt(0.75), s13 = std::sqrt(1./3.);
      double const center[2] = {0.5*edge, 0.5*s13*edge};
#endif

      for(int i = 0; i <= edge; ++i) { // must run serial due to non-local write operations
          for(int j = 0; j <= i; ++j) { // triangular loop
#ifdef GEOMETRIC
              int in, jn;
              {
                  // start from integer vector {i, j}
                  // apply base vectors to find position w.r.t. 00
                  double const v[2] = {i - 0.5*j, s34*j};
  #ifdef DEBUG
                  // check if the inversion works
                  double const ij_v[2] = {v[0] + s13*v[1], 2*s13*v[1]}; // this should match {i,j}
    //               printf("%g %g\n", ij_v[0] - i, ij_v[1] - j);
                  assert(std::abs(ij_v[0] - i) < 1e-13);
                  assert(std::abs(ij_v[1] - j) < 1e-13);
  #endif // debug
                  // subtract the center of the triangle
                  double const vc[2] = {v[0] - center[0], v[1] - center[1]};
                  // rotated by 120 degrees
                  double const wc[2] = {-.5*vc[0] + s34*vc[1], -s34*vc[0] - 0.5*vc[1]};
                  // add the center again
                  double const w[2] = {wc[0] + center[0], wc[1] + center[1]};
                  // convert into integer space
                  double const ij_w[2] = {w[0] + s13*w[1], 2*s13*w[1]};
                  // convert to integer values
                  in = std::round(ij_w[0]);
                  jn = std::round(ij_w[1]);
                  printf("# rotation (%d, %d) --> (%d, %d)\n", i,j, in,jn);
                  assert(0 <= in);
                  assert(in <= edge);
                  assert(0 <= jn);
                  assert(jn <= edge);
                  // confirm the integer formula below
                  assert(in == edge - i + j);
                  assert(jn == edge - i);
              }
#else
              // rotate by 120 degrees
              int const in = edge - i + j;
              int const jn = edge - i;
#endif // geometric
              assert(jn <= in);

#ifdef GEOMETRIC
              int im, jm;
              {
                  // start from integer vector {i, j}
                  // apply base vectors to find position w.r.t. 00
                  double const v[2] = {i - 0.5*j, s34*j};
                  // subtract the center of the triangle
                  double const vc[2] = {v[0] - center[0], v[1] - center[1]};
                  // mirror
                  double const wc[2] = {-vc[0], vc[1]};
                  // add the center again
                  double const w[2] = {wc[0] + center[0], wc[1] + center[1]};
                  // convert into integer space
                  double const ij_w[2] = {w[0] + s13*w[1], 2*s13*w[1]};
                  // convert to integer values
                  im = std::round(ij_w[0]);
                  jm = std::round(ij_w[1]);
                  printf("# mirror (%d, %d) --> (%d, %d)\n", i,j, im,jm);
                  assert(0 <= im);
                  assert(im <= edge);
                  assert(0 <= jm);
                  assert(jm <= edge);
                  assert(jm <= im);
                  // confirm the integer formula below
                  assert(im == edge - i + j);
                  assert(jm == j);
              }
#else
              // mirror
              int const im = edge - i + j;
              int const jm = j;
#endif // geometric
              assert(jm <= im);

              // 2*3 operation:
              //    1
              //    R
              //    M
              //    R*R
              //    R*M
              //    R*R*M

              int const i00 = i               , j00 = j;          // 1
              int const i0M = edge - i00 + j00, j0M = j00;        // M
              int const i10 = edge - i00 + j00, j10 = edge - i00; // R
              int const i1M = edge - i10 + j10, j1M = j10;        // R*M
              int const i20 = edge - i10 + j10, j20 = edge - i10; // R*R
              int const i2M = edge - i20 + j20, j2M = j20;        // R*R*M
              assert(j00 <= i00);
              assert(j10 <= i10);
              assert(j20 <= i20);
              assert(j0M <= i0M);
              assert(j1M <= i1M);
              assert(j2M <= i2M);

              // check that we do not mix levels
              auto const level_now = level(i,j);
              assert(level_now == level(i10,j10));
              assert(level_now == level(i20,j20));
              assert(level_now == level(i0M,j0M));
              assert(level_now == level(i1M,j10));
              assert(level_now == level(i2M,j2M));

              // now transfer the multiplicity
              auto & mr = multi(i,j); // add them up here
              { auto & m = multi(i10,j10); auto const mc{m}; m -= mc; mr += mc; }
              { auto & m = multi(i20,j20); auto const mc{m}; m -= mc; mr += mc; }
              { auto & m = multi(i0M,j0M); auto const mc{m}; m -= mc; mr += mc; }
              { auto & m = multi(i1M,j1M); auto const mc{m}; m -= mc; mr += mc; }
              { auto & m = multi(i2M,j2M); auto const mc{m}; m -= mc; mr += mc; }

          } // j
      } // i

      // mirror operation ee maps:
      //    00 <--> e0
      // implemented as (i,j) --> (e-i+j,j)

      // mirror (0, 0) --> (4, 0)
      // mirror (1, 0) --> (3, 0)
      // mirror (1, 1) --> (4, 1)
      // mirror (2, 0) --> (2, 0)
      // mirror (2, 1) --> (3, 1)
      // mirror (2, 2) --> (4, 2)
      // mirror (3, 0) --> (1, 0)
      // mirror (3, 1) --> (2, 1)
      // mirror (3, 2) --> (3, 2)
      // mirror (3, 3) --> (4, 3)
      // mirror (4, 0) --> (0, 0)
      // mirror (4, 1) --> (1, 1)
      // mirror (4, 2) --> (2, 2)
      // mirror (4, 3) --> (3, 3)
      // mirror (4, 4) --> (4, 4)


      // now generate all 2*3 operations as products

      std::vector<uint32_t> hm60(Level + 1, 0);
      std::vector<uint32_t> hist(128, 0);
      for(int i = 0; i <= edge; ++i) {
          for(int j = 0; j <= i; ++j) { // triangular loop
              ++hist[multi(i,j)];
              if (60 == multi(i,j)) ++hm60[level(i,j)];
          } // j
      } // i


      if (0) { // scope: eval
          if (echo > 0) printf("\n# %s: Histogramm of multiplicity 60 levels: (up to +Level=%d):\n", __func__, Level);
          for(size_t m = 0; m < hm60.size(); ++m) {
              if (hm60[m] > 0 && echo > 0) printf("# %d %d\n", m, hm60[m]); // we find hm60[L] == 2^{L - 1}
          } // m
      } // scope

      { // scope: eval
          if (echo > 0) printf("\n# %s: Histogramm of multiplicities (+Level=%d):\n", __func__, Level);
          size_t sum{0}, msum{0};
          for(size_t m = 0; m < hist.size(); ++m) {
              sum += hist[m];
              msum += m*hist[m];
              if (hist[m] > 0 && echo > 0) printf("# %d %d\n", m, hist[m]);
          } // m
          auto const reference = n_ico_vertices(Level);
          if (echo > 0) printf("# Histogramm checksum %ld reference %ld, sum %ld\n", msum, reference, sum);
          if (reference != msum) warn("total number of vertices %ld deviates from %ld", msum, reference);
          if (hist[12] != 1) warn("number of base points must be 10 + 2");
          if (Level > 0 && hist[30] != 1) warn("expect 30 edge-center points");
          if (Level > 1 && hist[60] < 1) warn("expect half-edge points with multiplicity 60");
          if (Level > 2 && hist[120] < 1) warn("number of face points must be 120");
          if (echo > 0) printf("# %s: L, m120, m60, n    %d %d %d %ld\n", __func__, Level, hist[120], hist[60], reference);
      } // scope

      return 0;
  } // test_triangle


  status_t inline test_base_point_distances(int const echo=0) {
      auto const SegmentDistance = 1./std::sin(Degrees72);
      if (echo > 7) printf("# %s: segment distance %.9f\n", __func__, SegmentDistance);
      double maxdev_norm{0}, maxdev_dist{0}, maxdev_pole{0};
      for(int i10 = 0; i10 < 10; ++i10) {
          double v[3]; base_point(v, i10);
          maxdev_norm = std::max(maxdev_norm, std::abs(1 - norm2(v)));
          for(int j10 = -2; j10 <= 2; ++j10) {
              if (0 != j10) {
                  double w[3]; base_point(w, i10 + j10);
                  add_product(w, 3, v, -1.); // subtract v from w
                  auto const dev = std::abs(pow2(SegmentDistance) - norm2(w));
                  maxdev_dist = std::max(maxdev_dist, dev);
              } // nonzero
          } // j10
          // measure distance from closest pole: north (south) pole for i10 even (odd)
          v[2] -= (1 - 2*(i10 % 2)); // subtract pole coordinates (only z coordinate)
          maxdev_pole = std::max(maxdev_pole, std::abs(pow2(SegmentDistance) - norm2(v)));
      } // i10
      if (echo > 5) printf("# %s: largest deviation from unit sphere is %.1e and from dodecahedron is %.1e and %.1e (poles)\n", 
                              __func__, std::sqrt(maxdev_norm), std::sqrt(maxdev_dist), std::sqrt(maxdev_pole));
      return (maxdev_norm > 1e-15);
  } // test_base_point_distances

  status_t inline test_generate(int const echo=1, unsigned const Level=4) {
      size_t const tL = rhomb_edge(Level);
      view4D<double> vtx(10, tL + 1, tL + 1, 4);
      return generate(vtx, Level, echo);
  } // test_generate

  status_t inline all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_base_point_distances(echo);
      unsigned const Level = control::get("Level", 6.);
      stat += test_generate(echo, Level);
//    for(unsigned L = 0; L <= Level; ++L)
      stat += test_triangle(echo, Level);
      return stat;
  } // all_tests

} // namespace icogrid
