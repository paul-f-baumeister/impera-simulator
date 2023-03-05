#pragma once

#include <cstdio> // std::printf, ::snprintf

#include "status.hxx" // status_t

namespace color {

  char constexpr esc = 27; //
  char const def[4] = {27,'[','m',0};

//   character(len=*), parameter, public :: def = esc//"[m" ! fall back to default font and background color
// 
//   ! standard font colors
//   character(len=*), parameter, public :: black = esc//"[30m", white = esc//"[37m", &
//     red = esc//"[31m", green = esc//"[32m", yellow = esc//"[33m", &
//     blue = esc//"[34m", magenta = esc//"[35m", cyan = esc//"[36m"
// 
//   ! standard background colors
//   character(len=*), parameter, public :: on_black = esc//"[40m", on_white = esc//"[47m", &
//     on_red = esc//"[41m", on_green = esc//"[42m", on_yellow = esc//"[43m", &
//     on_blue = esc//"[44m", on_magenta = esc//"[45m", on_cyan = esc//"[46m"
//   
//   character(len=*), parameter, public :: blink = esc//"[5m" ! font attribute blinking
// 
//   character(len=*), parameter, public :: green4 = esc//"[38;5;22m", & ! darker green ! esc//"[38;2;0;127;0m" ! KDE extension to 24-bit RGB
//                                       on_green4 = esc//"[48;5;22m"    ! darker green for the background
// 
//   character(len=*), parameter, public :: bold = esc//"[1m" ! font weight bold
//   character(len=*), parameter, public :: underline = esc//"[4m" ! singly underlined text
//   character(len=*), parameter, public :: invert = esc//"[7m" ! invert colors
// 
// !   character(len=*), parameter, public :: frame = esc//"[51m", encircle = esc//"[52m", esc//"[53m" ! do not work on Linux/OSX
//   
// 
// !   interface color
// !     module procedure color_font, color_background
// !   endinterface
// 
//   contains
// 
// !   character(len=16) function color_font(font, bg) result(str)
// !     character(len=*), intent(in)           :: font ! font color by name
// !     integer(kind=1),  intent(in), optional :: bg(3) ! background color in (0..5,0..5,0..5)
// !     selectcase( font )
// !     case( "white", "White" )    ; str = esc//"[37m"
// !     case( "black", "Black" )    ; str = esc//"[30m"
// !     case( "red", "Red" )        ; str = esc//"[31m"
// !     case( "green", "Green" )    ; str = esc//"[32m"
// !     case( "yellow", "Yellow" )  ; str = esc//"[33m"
// !     case( "blue", "Blue" )      ; str = esc//"[34m"
// !     case( "magenta", "Magenta" ); str = esc//"[35m"
// !     case( "cyan", "Cyan" )      ; str = esc//"[36m"
// !     case( "def", "default", "" ); str = esc//"[m"
// !     case default                ; str = esc//"[m"
// !     endselect
// !   endfunction
// ! 
// !   character(len=16) function color_background(font, bg) result(str)
// !     integer(kind=1),  intent(in)           :: font(3) ! font color in (0..5,0..5,0..5)
// !     character(len=*), intent(in), optional :: bg ! background color by name
// ! 
// !   endfunction
// 
// ! In 256 color mode (ESC[38;5;<fgcode>m and ESC[48;5;<bgcode>m), the color-codes are the following:[citation needed]
// ! 
// !  0x00-0x07:  standard colors (as in ESC [ 30..37 m)
// !  0x08-0x0f:  high intensity colors (as in ESC [ 90..97 m)
// !  0x10-0xe7:  6*6*6=216 colors: 16 + 36*r + 6*g + b (0<=r,g,b<=5)
// !  0xe8-0xff:  grayscale from black to white in 24 steps
//   
// 
// ! \e[31m --> red font color
// ! \e[31;1m --> red bold font color
// ! \e[41m --> red background color
// 
// ! Intensity       0       1       2       3       4       5       6       7
// ! Normal          Black   Red     Green   Yellow  Blue    Magenta Cyan    White
// ! Bright          Black   Red     Green   Yellow  Blue    Magenta Cyan    White
// 

  template<unsigned Nchars=8>
  class string_t {
  private:
      char c[Nchars];
  public:
      string_t() { if (Nchars) c[0] = '\0'; };
      string_t(char const *string) { std::snprintf(c, Nchars, string); }
      char const* c_str() const { return (Nchars) ? c : nullptr; }
      char*        data()       { return (Nchars) ? c : nullptr; }
  }; // class string_t

//   typedef string_t<8>  string7_t;
  typedef string_t<16> string15_t;

  inline string15_t colorchar(float const rd, float const gr, float const bl) { // rgb <= 1
      int const i = 16 + int(std::round(rd*5) + 6*(std::round(gr*5) + 6*std::round(bl*5)));
      string15_t s; std::snprintf(s.data(), 15, "%c[48;5;%im", esc, i);
      return s;
  } // colorchar

  inline string15_t colorchar(float const rgb[3]) { return colorchar(rgb[0], rgb[1], rgb[2]); }
  
  status_t inline test(int const echo=0) {
    
    string15_t c;

    for (int ib = 5; ib >= 0; --ib) {
      for (int ig = 5; ig >= 0; --ig) {
        for (int ir = 5; ir >= 0; --ir) {
            std::printf("%s ", colorchar(ir/5., ig/5., ib/5.).c_str());
        } // ir
      } // ig
      std::printf("%s\n", def); // reset and newline 
    } // ib
    std::printf("%s\n", def); // reset and newline 

    
//     for (int i = 0; i < 256; ++i) {
//         char c[16]; std::snprintf(c, 16, "%c[48;5;%im", esc, i);
//         char const x = ((i - 16) % 43) ? ' ' : 'x';
//       std::printf("(2A)",advance="no") trim(c),x
//       if(any([7,15,(15+k*36, k=0,5),231,255]==i)) std::printf("(2A)",advance="yes") c(1:2),"m" 
//     } // i
// 
//     do i = 16, 255, 1+6+36
//       std::printf(unit=c,fmt="(A,I0,A)") "e[48;5;",i,"m" ; c(1:1) = esc
//       std::printf("(2A)",advance="no") trim(c)," " 
// !       if(any([7,15,(15+k*36, k=0,5),231,255]==i)) std::printf("(2A)",advance="yes") c(1:2),"m" 
//     enddo

//     std::printf("%s%s%s\n", red,  "this should appear in red color", def);
//     std::printf("%s%s%s\n", green, "this should appear in green color", def);
//     std::printf("%s%s%s\n", blue, "this should appear in blue color", def);
//     std::printf("%s%s%s%s\n", blue, on_white, "this should appear in blue color on white background", def);
//     std::printf("%s%s%s%s%s\n", bold, blue, on_white, "this should appear in bold blue color on white background", def);
//     std::printf("%s%s%s%s%s\n", blink, red, on_white, "this should appear in blinking red color on white background", def);
//     std::printf("%s\n", "this should appear in default color");
//     std::printf("%s%s%s%s\n", green4, on_white, "this should appear in dark green color on white background", def);
//     std::printf("%s%s%s%s\n", green4, on_green, "      ", def);
//     std::printf("%s%s%s%s%s\n", green4, on_green, " T  @ ", def, " a tree on grass and some bushes");
//     std::printf("%s%s%s%s\n", green4, on_green, "   v  ", def);
//     std::printf("%s%s%s%s\n", green4, on_green, "v    @", def);
//     std::printf("%s%s%s\n", bold, "this should appear in bold default color", def);
//     std::printf("%s%s%s\n", underline, "this should appear underlined in default color", def);
//     std::printf("%s%s%s\n", invert, "this should appear in inverted default colors", def);
//     std::printf("%s%s%s\n", frame, "this should appear framed in default color", def);
//     std::printf("%s%s%s\n", encircle, "this should appear encircled in default color", def);
//     std::printf("%s%s%s\n", overline, "this should appear overlined in default color", def);

      return 0;
  } // test 

  
  inline
  status_t all_tests(int const echo=1) { return test(echo); }
  
} // namespace color
