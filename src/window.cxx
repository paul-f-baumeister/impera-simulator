#include <cstdio> // std::printf
#include <cassert> // assert
#include <algorithm> // std::min
#include <cmath> // M_PI

#include "window.hxx"

#include "data_view.hxx" // view4D<T>
#include "inline_tools.hxx" // set, add_product, intpow
#include "status.hxx" // status_t
#include "icogrid.hxx" // ::n_ico_vertices
#include "icomap.hxx" // ::GridCell_t
#include "control.hxx"  // ::get
#include "impera.hxx" // ::read_resource_file, ::run_one_day
// #include "impera.hxx" // ::Impera_t<real_t,Nspecies>

#ifdef HAS_WINDOW
// c++ -std=c++11 ... -framework GLUT -framework OpenGL -D GL_SILENCE_DEPRECATION

  // #define __MacOSX__
  #ifdef  __MacOSX__
    #include <OpenGL/gl.h>
    #include <OpenGL/glu.h>
    #include <OpenGL/glext.h>
    #include <GLUT/glut.h>
  #else  // __MacOSX__
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glext.h>
    #include <GL/glut.h>
  #endif // __MacOSX__

#endif // HAS_WINDOW

namespace window {

  // conversion factors for internal angle unit
  double constexpr i162rad = M_PI/32768.;
  double constexpr i162deg = 180./32768.;
  double constexpr deg2i16 = 32768./180.;
//   double constexpr rad2i16 = 32768./M_PI;
//   double constexpr rad2deg = 180./M_PI;
  double constexpr deg2rad = M_PI/180.;

  float transfer2float(uint32_t const i) { union{ uint32_t i; float f; } u{i}; return u.f; }
  uint32_t transfer2int32(float const f) { union{ float f; uint32_t i; } u{f}; return u.i; }

  // global quantities
  static int echo;
  static view4D<float> vertex0; // vertices of the unrotated wireframe model of the icosahedral sphere grid
  static view4D<float> vertex;  // vertices of the wireframe model of the icosahedral sphere grid
  static view3D<float> colors;  // actual colors of the triangles (view3D instead of view2D to query dim1)
  static int16_t itheta_now=9, iphi_now=9; // range of int16_t is [-32768, 32767]
  constexpr float zoom_max = 9;
  static float zoom = zoom_max;
  static float home_coords[2] = {99, 0}; // in degrees
  static GLfloat click_position[3][2][2]; // 2D positions where last clicked [left/middle/right][down/up][x/y]
  static int16_t   click_angles[3][2][2]; // 2D angles when last clicked [left/middle/right][down/up][theta/phi]
  static float     click_zoom  [3][2];
  static bool  mouse_button_is_down[3] = {false, false, false};
  static void* simulation = nullptr;
  static unsigned Nspecies = 0;

  template <typename real_t, int n=3>
  void matrix_matrix(real_t axb[n][n+1], real_t const a[n][n+1], real_t const b[n][n+1]) {
      for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
              real_t tmp{0};
              for (int k = 0; k < n; ++k) {
                  tmp += a[i][k] * b[k][j];
              } // k
              axb[i][j] = tmp;
          } // j
          axb[i][n] = 0;
      } // i
  } // matrix_matrix


  bool rotate(int itheta=0, int iphi=0, double psi=0) {
      itheta = std::min(std::max(-16383, itheta), 16383); // (-pi/2, pi/2), we do not allow to turn over the pole
          // we exclude the value 16384 which is pi/2 to avoid singularities as std::cos(16383*i162rad) =~= 1e-4
      if (itheta == itheta_now && iphi == iphi_now) return false;

      // convert internal int16_t representation of angles into radians
      double const theta = itheta*i162rad; // snap to grid
      double const phi = (-90 - 4)*deg2rad - iphi*i162rad; // 4 degrees from the map shift so that SouthAmerica fits into one rhomb

      // rotate vertex0 --> vertex
      auto const ca = std::cos(phi), // alpha
                 sa = std::sin(phi),
                 cb = std::cos(theta), // beta
                 sb = std::sin(theta);
//                  cg = std::cos(psi),   // gamma (roll angle?)
//                  sg = std::cos(psi);
//       double const matrix[3][4] = // general 3D rotation matrix with yaw, pitch and roll
//                  {{sa*sb*cg - ca*sg, ca*sb*cg + sa*sg, cb*cg,    0},
//                   {sa*sb*sg + ca*cg, ca*sb*sg - sa*cg, cb*sg,    0},
//                   {       sa*cb    ,        ca*cb    ,  -sb ,    0}};
      double matrix[3][4];
      {
          double const matrix_phi[3][4]   = {{ ca, -sa,  0 ,    0}, // turn around the z-axis
                                             { sa,  ca,  0 ,    0},
                                             { 0 ,  0 ,  1 ,    0}};

          double const matrix_theta[3][4] = {{ 1 ,  0 ,  0 ,    0}, // turn around the new x-axis
                                             { 0 ,  cb, -sb,    0},
                                             { 0 ,  sb,  cb,    0}};
          double matrix_theta_phi[3][4];
          matrix_matrix(matrix_theta_phi, matrix_theta, matrix_phi);
          // switch y <--> -z as in OpenGL, +y is upward
          double const matrix_swap[3][4] = {{ 1,  0,  0,    0},
                                            { 0,  0,  1,    0},
                                            { 0, -1,  0,    0}};
          matrix_matrix(matrix, matrix_swap, matrix_theta_phi);
      }
      // also 
      int const tL = vertex0.dim1() - 1;
      for (int i10 = 0; i10 < 10; ++i10) {
          for (int iSE = 0; iSE <= tL; ++iSE) { // south east direction
              for (int iNE = 0; iNE <= tL; ++iNE) { // north east direction
                  auto const *const vtx0 = vertex0(i10,iSE,iNE); // unrotated vectors
                  for (int i = 0; i < 3; ++i) {
                      double tmp{0};
                      for (int j = 0; j < 3; ++j) {
                          tmp += matrix[i][j] * vtx0[j];
                      } // j
                      vertex(i10,iSE,iNE,i) = tmp;
                  } // i
                  vertex(i10,iSE,iNE,3) = vertex0(i10,iSE,iNE,3); // warning, may not be regular floating point format
              } // iNE
          } // iSE
      } // i10

      itheta_now = itheta; iphi_now = iphi; // set global variables

      if (echo) std::printf("# %s (theta, phi) = (%c %.3f, %c %.3f) degrees, (%d, %d)/2^16 and (%9.6f, %9.6f) rad\n",
                          __func__, (itheta_now < 0)?'S':'N', std::abs(itheta_now)*i162deg, 
                                    (iphi_now   < 0)?'W':'E', std::abs(iphi_now  )*i162deg,
                                    itheta_now, iphi_now, itheta_now*i162rad, iphi_now*i162rad);
      return true;
  } // rotate


  inline float zoom_distance(int const z) { return std::pow(1.125f, z + 1); } // in units of the earth radius --> at least 825 km over ground

  bool change_zoom(float const new_zoom) {
      auto const old_zoom = zoom;
      zoom = std::min(std::max(0.f, new_zoom), zoom_max);
      bool const changed = (old_zoom != zoom);
      if (echo) std::printf("# %s%s from %.3f (%.6f) to %.3f (%.6f)\n", changed?"":"cannot ",__func__,
          zoom_distance(old_zoom), old_zoom, zoom_distance(zoom), zoom);
//       if (changed) {
//           GLfloat model[16];
//           set(model, 16, GLfloat(0));
//           glGetFloatv(GL_MODELVIEW_MATRIX, model);
//           std::printf("# GL_MODELVIEW_MATRIX:"); for (int i = 0; i < 16; ++i) std::printf("%s%12.6f", (i&3)?" ":"\n#", model[i]); std::printf("\n");
//       }
      return changed;
  } // change_zoom

  int16_t angle_of_zoom(float const z) { return std::pow(2.f, 4 + z); }

  status_t create_wireframe(unsigned const Level=3, int const echo=11) {
      window::echo = echo; // copy verbosity level
      status_t stat(0);
      { // scope: generate the vertices
          int const tL = icogrid::rhomb_edge(Level);
          view4D<double> vtx(10, tL + 1, tL + 1, 4);
          stat += icogrid::generate(vtx, Level, echo);
          if (0 != stat) {
              warn("failed to generate icogrid for Level=%d status=%i", Level, int(stat));
              return stat;
          }
          vertex0 = view4D<float>(10, tL + 1, tL + 1, 4, 0.f);
          for (int i10 = 0; i10 < 10; ++i10) {
              for (int iSE = 0; iSE <= tL; ++iSE) { // south east direction
                  for (int iNE = 0; iNE <= tL; ++iNE) { // north east direction
                      for (int d = 0; d < 3; ++d) {
                          vertex0(i10,iSE,iNE,d) = vtx(i10,iSE,iNE,d);
                      } // d
                      // store the ico_index in the 4-th component
                      auto const index = icogrid::ico_index_wrap(Level, i10, iSE, iNE);
                      uint32_t const index32 = index;
                      assert(index == index32);
                      vertex0(i10,iSE,iNE,3) = transfer2float(index32); // warning, may not be regular floating point format
                  } // iNE
              } // iSE
          } // i10

          colors = view3D<float>(1, icogrid::n_ico_vertices(Level), 3, .25f); // initialize color buffer all triangles dark grey
          vertex = view4D<float>(10, tL + 1, tL + 1, 4, 0.f);
          rotate(); // transfer vertex0 into vertex, initialize (theta_now,phi_now)
          
      } // scope
      return stat;
  } // create_wireframe

  
  void where_was_the_click(int const x, int const y, GLfloat xy[2]=nullptr) {
      auto const width  = glutGet(GLUT_WINDOW_WIDTH); // ToDo: move these into global variables that get updated on resize
      auto const height = glutGet(GLUT_WINDOW_HEIGHT);
      if (xy) {
          xy[0] = x/(0.5*width) - 1.;
          xy[1] = 1. - y/(0.5*height);
      } else {
          auto const min_wh = std::max(1, std::min(width, height));
          auto const xr = (x -  width*.5)/min_wh;
          auto const yr = (y - height*.5)/min_wh;
          std::printf("# click at (x=%d, y=%d) --> (%g, %g)\n", x, y, xr, yr);
      }
      // at distance 3.247 earth radii, the earth circumference is at xr == yr == +/-0.39
      // at distance 2.887 earth radii, the earth circumference is at xr == yr == +/-0.446
      // at distance 2.027 earth radii, the earth circumference is at xr == yr == +/-0.6845
  } // where_was_the_click
  
  void mouse_click(int const button, int const state, int const x, int const y) {
      // int constexpr WHEEL_UP = 3, WHEEL_DOWN = 4; // https://stackoverflow.com/questions/14378/using-the-mouse-scrollwheel-in-glut seems not to work on OSX
      int const button012 = (GLUT_LEFT_BUTTON   == button) ? 0 : (
                            (GLUT_MIDDLE_BUTTON == button) ? 1 : (
                            (GLUT_RIGHT_BUTTON  == button) ? 2 : 3));
      char const botton_string[][8] = {"left", "middle", "right", "?"};
      int const state01  = (GLUT_DOWN == state) ? 0 : ((GLUT_UP == state) ? 1 : 2);
      int const down0up1 = (GLUT_DOWN == state) ? 0 : 1;
      char const state_string[][8] = {"press", "release", "?"};
      std::printf("# Mouse %s %s at x=%d y=%d\n", botton_string[button012], state_string[state01], x, y);
      // beware: x,y can be negative when GLUT_UP == state but not when GLUT_DOWN == state
      if (GLUT_DOWN == state) assert(x >= 0 && y >= 0);
      where_was_the_click(x, y);
      if (button012 < 3) {
          // store the current angles and zoom
          click_angles[button012][down0up1][0] = itheta_now;
          click_angles[button012][down0up1][1] = iphi_now;
          click_zoom[button012][down0up1] = zoom;
          where_was_the_click(x, y, click_position[button012][down0up1]);
          if ((0 == down0up1) && mouse_button_is_down[button012]) {
              warn("%s mouse button was already down", botton_string[button012]);
          }
          mouse_button_is_down[button012] = (GLUT_DOWN == state);
      }
      glutPostRedisplay(); // redraw
  } // mouse_click

  void mouse_motion(int const x, int const y) {
      // std::printf("# Mouse motion x=%d y=%d\n", x, y); // tons of output when you move the mouse, not only when pressed
      bool changed{false};
      if (mouse_button_is_down[2]) {
          int const drag_button = 2;
          // free earth rotation mode
          GLfloat xy[2]; where_was_the_click(x, y, xy);
          auto const f = 10*deg2i16;
          int16_t const itheta_new = click_angles[drag_button][0][0] + f*(click_position[drag_button][0][1] - xy[1]);
          int16_t const iphi_new   = click_angles[drag_button][0][1] + f*(click_position[drag_button][0][0] - xy[0]);
          std::printf("# Mouse drag at x=%d y=%d --> new angles (%g, %g)\n", x, y, itheta_new*i162deg, iphi_new*i162deg);
          changed = rotate(itheta_new, iphi_new);
      } else if (mouse_button_is_down[1]) { // middle button
          int const drag_button = 2;
          // free earth rotation mode
          GLfloat xy[2]; where_was_the_click(x, y, xy);
          auto const f = 1;
          auto const new_zoom = click_zoom[drag_button][0] + f*(click_position[drag_button][0][1] - xy[1]);
          std::printf("# Mouse drag at y=%d --> new zoom %g\n", y, new_zoom);
          changed = change_zoom(new_zoom);
      }
      
      // else { std::printf("# Mouse move at x=%d y=%d\n", x, y); }
      if (changed) glutPostRedisplay(); // redraw
  } // mouse

  inline bool focus_on(double const theta_degree, double const phi_degree, float const zoom_level=zoom_max) {
      return change_zoom(zoom_level) || rotate(theta_degree*deg2i16, phi_degree*deg2i16);
  } // focus_on
  
  static void key_stroke(int const key, int const x, int const y) {
      // beware: x,y can be negative
      auto const angle = angle_of_zoom(zoom);
      bool changed{false};
      switch (key) {
          case GLUT_KEY_F1:
          case GLUT_KEY_F2:
          case GLUT_KEY_F3:
          case GLUT_KEY_F4:
          case GLUT_KEY_F5:
          case GLUT_KEY_F6:
          case GLUT_KEY_F7:
          case GLUT_KEY_F8: // also the delete key on OSX
          case GLUT_KEY_F9:
          case GLUT_KEY_F10:
//        case GLUT_KEY_F11: // cannot be used on OSX
          case GLUT_KEY_F12:
              if (echo > 9) std::printf("# %s key=F%d x=%d y=%d\n", __func__, key, x, y);
          break;

          case GLUT_KEY_LEFT      : changed = rotate(itheta_now, iphi_now - angle); break;
          case GLUT_KEY_RIGHT     : changed = rotate(itheta_now, iphi_now + angle); break;
          case GLUT_KEY_UP        : changed = rotate(itheta_now + angle, iphi_now); break;
          case GLUT_KEY_DOWN      : changed = rotate(itheta_now - angle, iphi_now); break;

          case GLUT_KEY_PAGE_UP   : changed = change_zoom(int(zoom - 1)); break;
          case GLUT_KEY_PAGE_DOWN : changed = change_zoom(int(zoom + 1)); break;

          case GLUT_KEY_END       :
          case GLUT_KEY_INSERT    :
              if (echo > 9) std::printf("# %s key=GLUT_KEY#%d x=%d y=%d\n", __func__, key, x, y);
          break;

          // focus on special regions of the world
          case 'A': changed = focus_on(30, 120); break; // center on Asia
          case 'I': changed = focus_on(20, 78, 4); break; // center on India
          case 'E': changed = focus_on(48, 16); break; // center on Europe
          case 'N': changed = focus_on(40, -100); break; // center on NorthAmerica
          case 'S': changed = focus_on(-13, -58); break; // center on SouthAmerica
          case 'O': changed = focus_on(-25, 136); break; // center on Oceania
          // Cities
          case 'a': changed = focus_on(52.5, 13.4, 0); break; // center on Berlin
          case 'b': changed = focus_on(52.5, 13.4, 0); break; // center on Berlin
          case 'c': changed = focus_on(30, 31.2, 0); break; // center on Cairo
//           case 'd': changed = focus_on(28.6, 77.2, 0); break; // center on newDelhi
//           case 'h': changed = focus_on(21.3, -157.9, 0); break; // center on Honolulu/Hawaii
//           case 'j': changed = focus_on(-6.2, 106.8, 0); break; // center on Jakarta
//           case 'k': changed = focus_on(24.9, 67, 0); break; // center on Karachi
          case 'M': changed = focus_on(14.6, 121, 0); break; // center on Manila
          case 'm': changed = focus_on(40, -4, 0); break; // center on Madrid
          case 'n': changed = focus_on(41, -74, 0); break; // center on NewYork
          case 's': changed = focus_on(37.5, 127, 0); break; // center on Seoul, missing Shanghai
          case 't': changed = focus_on(35.7, 139.8, 0); break; // center on Tokyo

          case ' ': // space for start/stop
              if (echo > 17) std::printf("# %s key=[space] x=%d y=%d\n", __func__, x, y);
              if (simulation) {
                  // TODO spawn a thread and run while(startstop)
                  impera::run_one_day(Nspecies);
                  changed = true; // --> redraw
              } else warn("simluation has not been initialized");
          break;

          case GLUT_KEY_HOME:
              if (99 == home_coords[0]) {
                  home_coords[0] = 45;
                  warn("Home coordinates have not been set, use (%g, %g)", home_coords[0], home_coords[1]);
              }
              changed = focus_on(home_coords[0], home_coords[1], 5);
          break;

          case 127: // DELETE
              if (echo > 9) std::printf("# %s key=DEL x=%d y=%d\n", __func__, x, y);
          break;
          
          default:
              if (key < 32) {
                  if (echo > 9) std::printf("# %s key=ASCII#%d x=%d y=%d\n", __func__, key, x, y); // non printable characters
              } else if (key < 127) {
                  if (echo > 9) std::printf("# %s key=\'%c\' x=%d y=%d\n", __func__, char(key), x, y);
              }
          break;
      } // switch key
      if (changed) {
          glutPostRedisplay(); // redraw
      }
  } // key_stroke

  void points_to_go(); // forward declaration

  void timer(int const val) {
      glutPostRedisplay(); // redraw
      points_to_go();
  } // timer

  void points_to_go() {
      glutTimerFunc(1, timer, 0);
  } // points_to_go

  void show_vertex(float const pos[4], float const rgb[3]=nullptr) {
      if (!rgb) {
          bool constexpr colorful_sphere = false;
          if (colorful_sphere) {
              rgb = pos;
          } else {
              auto const ic = transfer2int32(pos[3]); // color index
              // std::printf("# color index = %d\n", ic);
              rgb = colors(0,ic);
          }
      }
      glColor3f(rgb[0], rgb[1], rgb[2]); // set the color
      glVertex3f(pos[0], pos[1], pos[2]); // place a vertex
  } // show_vertex

  void display(void) {

//           float const red  [] = {1.f, 0.f, 0.f};
//           float const green[] = {0.f, 1.f, 0.f}; // base colors
//           float const blue [] = {0.f, 0.f, 1.f};

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
      glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix
    
      // Render a sphere with a icosahredral grid
      glLoadIdentity();               // Reset the model-view matrix
      glTranslatef(0.f, 0.f, -zoom_distance(zoom)); // Move right and into the screen
      glBegin(GL_TRIANGLES);          // Begin drawing triangles

      int const tL = vertex.dim1() - 1;
      for (int i10 = 0; i10 < 10; ++i10) { // rhombs
          for (int iSE = 0; iSE < tL; ++iSE) { // south east direction
              for (int iNE = 0; iNE < tL; ++iNE) { // north east direction
                  // northern triangle
                  show_vertex(vertex(i10,iSE    ,iNE    ));
                  show_vertex(vertex(i10,iSE + 1,iNE + 1));
                  show_vertex(vertex(i10,iSE    ,iNE + 1));
                  // southern triangle
                  show_vertex(vertex(i10,iSE    ,iNE    ));
                  show_vertex(vertex(i10,iSE + 1,iNE    )); // e.g. the south pole is blue
                  show_vertex(vertex(i10,iSE + 1,iNE + 1));
              } // iNE
          } // iSE
      } // i10

      glEnd();                        // End of drawing triangles

      // draw the cursor
      if (mouse_button_is_down[0]) { 
          glBegin(GL_POINTS); // draw a single point
          glColor3f(1.f, 1.f, 1.f); // set the color pure white
          glPointSize(1);
//        std::printf("# Draw a point at x=%g y=%g\n", click_position[0][0][0], click_position[0][0][1]);
          glVertex2fv(click_position[0][0]); // normalized device coordinates
          glEnd();                        // End of drawing a point
          // point is drawn but the coordinate transformation seem unclear
      } // draw a point when clicked
      
      glutSwapBuffers(); // Swap the front and back frame buffers (double buffering)
    
  } // display


  void display3D(int const Level, float const *rgb) {
      if (!rgb)      return;
      if (Level < 0) return;
      auto const n = icogrid::n_ico_vertices(Level);
      if (colors.dim1() < n) {
          // re-allocate global array 'colors'
          colors = view3D<float>(1, n, 3, 0.f); // (view3D instead of view2D to query dim1)
      }
      set(colors.data(), n*3, rgb); // deep copy
      display();
  } // display3D


  /* Handler for window re-size event. Called back when the window first appears and
    whenever the window is re-sized with its new width and height */
  void reshape_window(GLsizei const width, GLsizei const height) {  // GLsizei for non-negative integer
      // Compute aspect ratio of the new window
      GLfloat const aspect = width/GLfloat(std::max(GLsizei(1), height));

      // Set the viewport to cover the new window
      glViewport(0, 0, width, height);

      // Set the aspect ratio of the clipping volume to match the viewport
      glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
      glLoadIdentity();             // Reset
      // Enable perspective projection with fovy, aspect, zNear and zFar
      gluPerspective(45.f, aspect, 0.1f, 100.f);
  } // reshape_window

  template<typename real_t>
  status_t set_color_map(view3D<real_t> const & inp, int const echo=9) {
      status_t stat(0);
      int const nr = colors.dim1();
      if (inp.dim1() == nr && inp.stride() >= 3) {
          real_t max_resources{0};
          for (int ir = 0; ir < nr; ++ir) {
              for (int rgb = 0; rgb < 3; ++rgb) {
                  auto const value = inp(0,ir,rgb);
                  auto const r = std::sqrt(std::sqrt(std::abs(value)));
                  colors(0,ir,rgb) = r;
                  max_resources = std::max(max_resources, r);
              } // rgb
          } // ir
          float const by_max = 1./max_resources;
          for (int ir = 0; ir < nr; ++ir) {
              for (int rgb = 0; rgb < 3; ++rgb) {
                  colors(0,ir,rgb) *= by_max; // scale
              } // rgb
          } // ir
      } else warn("dimension 1 do not agree, colors(1,%d,3) vs inp(?,%d,%d)", nr, inp.dim1(), inp.stride());
      return stat;
  } // set_color_map


  status_t get_resource_map(unsigned const Level=3, int const echo=9) {
      status_t stat(0);
      auto const ResourceFileBaseName  = control::get("ResourceFileBaseName", "mapfile_L");
      auto const ResourceFileExtension = control::get("ResourceFileExtension", ".dat");
      auto const ResourceFilePath      = control::get("ResourceFilePath", "../icogrid");
      char infilename[96];
      std::snprintf(infilename, 95, "%s/%s%i%s", ResourceFilePath, ResourceFileBaseName, Level, ResourceFileExtension);
      if (echo > 0) std::printf("# read resource map from '%s'\n", infilename);
      std::vector<GridCell_t> resources;
      stat += impera::read_resource_file(resources, infilename, echo);
//       int const nr = std::min(, colors.dim1());
//       float max_resources{0};
//       for (int ir = 0; ir < nr; ++ir) {
//           auto const r = std::sqrt(std::sqrt(resources[ir].resources));
//           colors(0, ir, 0) = r;
//           max_resources = std::max(max_resources, r);
//       } // ir
//       float const by_max = 1./max_resources;
//       for (int ir = 0; ir < nr; ++ir) {
//           colors(0, ir, 0) *= by_max; // scale
//       } // ir
      view3D<float> resource_map(1, resources.size(), 3, 0.f);
      for (size_t ir = 0; ir < resources.size(); ++ir) {
          int constexpr RED = 0;
          resource_map(0, ir, RED) = resources[ir].resources;
      } // ir
      stat += set_color_map(resource_map, echo);
      return stat;
  } // get_resource_map

  void run() {
      if (echo > 0) std::printf("# start run\n");
      glutMainLoop();
      if (echo > 0) std::printf("# run done\n");
  } // run
  
  status_t init(int argc, char *argv[]) {
      glutInit(&argc, argv);
      glutInitDisplayMode(GLUT_DOUBLE); // enable double buffered mode
      glutInitWindowSize(1024, 768); // initial window size
      glutInitWindowPosition(10, 10); // position of the window from the top left corner of the screen
      glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
      glutCreateWindow("Impera");
      glutDisplayFunc(display);
      glutSpecialFunc(key_stroke);
      glutReshapeFunc(reshape_window);       // Register callback handler for window re-size event
      glutMouseFunc(mouse_click);
//       glutPassiveMotionFunc(mouse_motion); // use __Passive__ if no button is pressed
      glutMotionFunc(mouse_motion);
      glClearColor(0, 0, 0, 1); // set background white

      glClearColor(0.f, 0.f, 0.f, 1.f); // Set background color to black and opaque
      glClearDepth(1.f);                   // Set background depth to farthest
      glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
      glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
      glShadeModel(GL_SMOOTH);   // Enable smooth shading
      glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
      glEnable(GL_CULL_FACE);  // show only front faces
      return 0;
  } // init

  status_t finalize() {
      // free dyamic memory in global quantities
      vertex0 = view4D<float>();
      vertex  = view4D<float>();
      colors  = view3D<float>();
      return 0;
  } // finalize

  status_t inline test_world_map(int const echo=1) {
      status_t stat(0);
      int const Level = control::get("Level", 3.); // icosahedral grid level
      stat += create_wireframe(Level, echo);
      stat += get_resource_map(Level, echo);
      stat += init(0, nullptr);
      control::set("RenderTime", "1");
      control::set("DisplayTime", "-1");
      Nspecies = control::get("Nspecies", 3.);
      simulation = impera::run_one_day(Nspecies);
      run();
      return stat;
  } // test_world_map

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_world_map(echo);
      stat += finalize();
      return stat;
  } // all_tests

} // namespace window
