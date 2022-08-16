#include <cstdio> // std::printf, ::snprintf
#include <cassert> // assert
#include <string> // std::string
#include <cstdlib> // std::atof
#include <map> // std::map<T1,T2>
// #include <utility> // std::pair<T1,T2>
#include <tuple> // std::tuple<T1,...,Tn>, ::get
#include <cstring> // std::strchr, ::strncpy
#include <cmath> // std::sqrt
#include <fstream> // std::fstream

#include "control.hxx" // default_echo_level
// declarations: set, get, command_line_interface, read_control_file, all_tests

#include "warnings.hxx" // warn

namespace control {

  int constexpr MaxNameLength = 64; // max variable name length

  inline void double2string(char buffer[32], double const number) {
      std::snprintf(buffer, 31, "%.16e", number);
  } // double2string

  inline double string2double(char const *string) {
      return std::atof(string);
  } // string2double


  int constexpr default_value_tag = 2e9;
  // hidden function: _manage_variables(echo, name, value) --> set
  //                  _manage_variables(echo, name, value, linenumber) --> set to default
  //                  _manage_variables(echo, name) --> get
  //                  _manage_variables(echo) --> show_variables
  char const* _manage_variables(
        int const echo
      , char const *name=nullptr
      , char const *value=nullptr
      , int const linenumber=0
  ) {

      static std::map<std::string, std::tuple<std::string,uint32_t,int32_t>> _map; // hidden archive

      if (name) {
          assert(nullptr == std::strchr(name, '=')); // make sure that there is no '=' sign in the name

          std::string const varname(name);
          if (value) {

              // set
              auto const newvalue = std::string(value);
              bool constexpr warn_redefines = true;
              if (warn_redefines) {
                  auto const & oldvalue = std::get<0>(_map[varname]);
                  auto const old = oldvalue.c_str();
                  bool const redefined = (old && '\0' != *old);
                  if (echo > 7) {
                      std::printf("# control sets \"%s\"", name);
                      if (redefined) std::printf(" from \"%s\"", old);
                      std::printf(" to \"%s\"\n", value);
                  } // echo
                  if (redefined) {
                      warn("variable \"%s\" was redefined from \"%s\" to \"%s\"", name, old, value);
                  } // redefined
              } // warn_redefines
              auto & tuple = _map[varname];
              std::get<0>(tuple) = newvalue;
              std::get<1>(tuple) = (default_value_tag == linenumber); // counter how many times this variable was evaluated
              std::get<2>(tuple) = linenumber; // store line number in input file 
                                               // or (if negative) command line argument number
              return value;

          } else { // value

              // get
              auto & tuple = _map[varname];
              auto const oldvalue = std::get<0>(tuple).c_str();
              ++std::get<1>(tuple); // increment reading counter
              if (echo > 7) std::printf("# control found \"%s\" = \"%s\"\n", name, oldvalue);
              return oldvalue;

          } // value

      } else { // name

          // show_variables
          if (echo) {
              int const show = echo;
              std::printf("\n# control.show=%d (1:minimal, 2:unused, 4:defaults, negative for details)\n", show);
              bool const show_unused  = std::abs(show) & 0x2; // all or only accessed ones
              bool const show_default = std::abs(show) & 0x4; // list also variables that are at their default value
              bool const show_details = (show < 0); // show access count and is_default
              std::printf("# control has the following variables defined:\n#\n");
              int listed{0};
              for (auto const & pair : _map) {
                  auto const times_used = std::get<1>(pair.second); // how many times was this value used?
                  if (show_unused || times_used > 0) {
                      int const line = std::get<2>(pair.second);
                      bool const is_default = (default_value_tag == line);
                      if (show_default || !is_default) {
                          auto const string = std::get<0>(pair.second);
                          double const numeric = string2double(string.c_str());
                          char buffer[32]; double2string(buffer, numeric);
                          if (show_details) {
                              std::printf("# used %dx, %s %d\t\t", times_used,
                                  is_default?"def ":((line > 0)?"argv":(line?"line":"set ")),
                                  is_default?0:std::abs(line));
                          }
                          if (string == buffer) {
                              // can be parsed as double, print with %g format
                              std::printf("# %s=%g\n", pair.first.c_str(), numeric);
                          } else {
                              std::printf("# %s=%s\n", pair.first.c_str(), string.c_str());
                          }
                          ++listed;
                      } // variable is at its default value
                  } // variable has been used or we want to show all variables defined
              } // pair
              std::printf("#\n# %d variables listed for control.show=%d\n", listed, show);
          } // show
          return nullptr;

      } // name

  } // _manage_variables

  void set(char const *name, char const *value, int const echo) {
      if (echo > 5) std::printf("# control::set(\"%s\", \"%s\")\n", name, value);
      assert(nullptr != value);
      _manage_variables(echo, name, value); // set
  } // set<string>

  char const* get(char const *name, char const *default_value) {
      int const echo = default_echo_level;
      auto const value = _manage_variables(echo, name); // query
      if (nullptr != value && '\0' != *value) {
          if (echo > 5) std::printf("# control::get(\"%s\", default=\"%s\") = \"%s\"\n", name, default_value, value);
          return value;
      } else {
          if (echo > 5) std::printf("# control::get(\"%s\") defaults to \"%s\"\n", name, default_value);
          return _manage_variables(echo, name, default_value, default_value_tag); // set to default
      }
  } // get<string>

  status_t show_variables(int const echo) {
      _manage_variables(echo);
      return 0;
  } // show_variables

  inline char* find_equal_sign(char const *string) { 
      return (char*)std::strchr(string, '=');
  } // find_equal_sign

  status_t command_line_interface(char const *statement, int const iarg, int const echo) {
      auto const equal = find_equal_sign(statement);
      if (nullptr == equal) {
          warn("ignored statement \"%s\", maybe missing \'=\'", statement);
          return 1; // error, no '=' sign given
      } else {
          auto const equal_char = equal - statement;
          assert('=' == statement[equal_char]);
          char name[MaxNameLength]; // get a mutable string
          std::strncpy(name, statement, std::min(MaxNameLength, int(equal_char))); // copy the statement up to '='
          name[equal_char] = '\0'; // delete the '=' sign to mark the name
          char const *value = equal + 1; // everything after the '=' marker
          if (echo > 7) std::printf("# control::set(statement=\"%s\") found name=\"%s\", value=\"%s\"\n", statement, name, value);
          _manage_variables(echo, name, value, iarg); // set the variable
          return 0;
      }
  } // command_line_interface

  void set(char const *name, double const value, int const echo) {
      /* This version of set is only used in the tests below */
      char buffer[32]; double2string(buffer, value);
      return set(name, buffer, echo);
  } // set<double>

  double get(char const *name, double const default_value) {
      char buffer[32]; double2string(buffer, default_value);
      return string2double(get(name, buffer));
  } // get<double>


  std::string left_trim(std::string const & s)  {
      std::string const WhiteSpaceChars = " \n\r\t\f\v";
      size_t const start = s.find_first_not_of(WhiteSpaceChars);
      return (start == std::string::npos) ? "" : s.substr(start);
  } // left_trim

//   std::string trim(std::string const & s)  {
//       size_t const end = s.find_last_not_of(WhiteSpaceChars);
//       return (end == std::string::npos) ? "" : s.substr(0, end + 1);
//   } // right_trim

  status_t read_control_file(char const *filename, int const echo) {
      status_t stat(0);
      char const CommentChar = '#'; // commented lines in control files
      char const EchoComment = '!'; // comments that should appear in the log

      if (nullptr == filename) {
          if (echo > 1) std::printf("# no control file passed\n");
          return stat; // 0
      }
      
      if ('\0' == *filename) {
          if (echo > 1) std::printf("# no control file given\n");
          return stat; // 0
      }

      std::ifstream infile(filename, std::ifstream::in);
      if (infile.fail()) {
          warn("Unable to open file '%s' for reading controls", filename);
          return -1;
      } // failed

      if (echo > 1) std::printf("\n# reading '%s' ...\n\n", filename);
      int linenumber{0}, ncomments{0}, nempty{0};
      std::string line;
      while (std::getline(infile, line)) {
          ++linenumber;
          if (echo > 18) std::printf("# %s:%d\t  %s\n", filename, linenumber, line.c_str());
          auto const tlin = left_trim(line);
          if (CommentChar == tlin[0]) {
              ++ncomments;
              if (echo > 9) std::printf("# %s:%d\t comment: %s\n",
                                           filename, linenumber, tlin.c_str());
              if (EchoComment == tlin[1]) {
                  if (echo > 0) std::printf("%s\n", tlin.c_str());
              }
          } else if ("" == tlin) {
              ++nempty;
              if (echo > 11) std::printf("# %s:%d\t is empty\n", filename, linenumber);
          } else {
              if (echo > 8) std::printf("# %s:%d\t  %s\n", filename, linenumber, tlin.c_str());
              auto const line_stat = command_line_interface(tlin.c_str(), -linenumber);
              if (line_stat) {
                  warn("failure parsing %s:%d \'%s\'", filename, linenumber, line.c_str());
              } else {
                 if (echo > 0) std::printf("# %s\n", tlin.c_str()); // show the valid commands
              } 
              stat += line_stat;
          }
      } // parse file line by line

      if (echo > 3) std::printf("# %s found %d comments in file '%s', status=%i\n\n",
                                   __func__, ncomments, filename, int(stat)); 

      return stat;
  } // read_control_file

  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_control(int const echo=9) {
      if (echo > 1) std::printf("\n# %s: %s\n", __FILE__, __func__);
      status_t stat(0);

      set("a", "5", echo); // use the string-setter routine
      if (echo > 1) std::printf("# a = %s\n", get("a", ""));

      set("b", 6., echo); // use the double-setter routine
      if (echo > 1) std::printf("# b = %s\n", get("b", ""));

      // or use the command line interface
      stat += command_line_interface("a=6", 99, echo); // launches warning about redefining "a"

      auto const a = get("a", "defaultA");
      if (echo > 1) std::printf("# a = %s\n", a);

      auto const b = get("b", "defaultB");
      if (echo > 1) std::printf("# b = %s\n", b);

      auto const c = get("c", "3.14");
      if (echo > 1) std::printf("# c = %s\n", c);

      auto const c_double = get("c", 3.1415);
      if (echo > 1) std::printf("# c<double> = %g\n", c_double);

      return stat;
  } // test_control

  status_t test_precision(int const echo=3, int const nmax=106) {
      // check if there are rounding errors arising from the 
      //    ASCII representation of double precision numbers
      if (echo > 2) std::printf("\n# %s: %s\n", __FILE__, __func__);
      status_t stat(0);
      double d{0.2}; // choose 1/5 which cannot be represented exactly in binary
      for (int i = 0; i < nmax; ++i) {
          set("d", d, echo); // warning about redefining "d" in the second iteration
          double const g = get("d", 1.);
          stat += (g != d);
          d *= std::sqrt(33/32.); // some irrational number close to 1
      } // i
      if (echo > 1) std::printf("# %s: for %i of %i cases double precision numbers are not retrieved\n", __func__, stat, nmax);
      return stat;
  } // test_precision

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_control(echo);
      stat += test_precision(echo*0);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace control
