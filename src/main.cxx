#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::abs
#include <vector> // std::vector
#include <string> // std::string

#include "simple_timer.hxx" // SimpleTimer
#include "warnings.hxx" // warn, ::show_warnings, ::clear_warnings
#include "control.hxx" // ::command_line_interface

#include "warnings.hxx" // ::all_tests
#include "inline_tools.hxx" // ::all_tests
#include "data_view.hxx" // ::all_tests
#include "control.hxx" // ::all_tests
#include "icomap.hxx" // ::all_tests
#include "impera.hxx" // ::all_tests
#include "bitmap.hxx" // ::all_tests
#include "icogrid.hxx" // ::all_tests
#include "window.hxx" // ::all_tests
#include "color.hxx" // ::all_tests

namespace main_tools {

  int run_unit_tests(char const *module=nullptr, int const echo=0) {
      bool const empty = (nullptr == module);
      std::string const module_name(empty ? "" : module);
      bool const show = ('?' == module_name[0]);
      bool const all = empty || show;
      if (echo > 0) {
          if (show) { printf("\n# show available module tests:\n"); } else
          if (all)  { printf("\n# run all tests!\n\n"); }
          else      { printf("# run unit tests for module '%s'\n\n", module_name.c_str()); }
      } // echo

      std::vector<std::pair<char const*,status_t>> results;
      { // testing scope
#define   add_module_test(NAME) \
          if (all || (0 == module_name.compare(#NAME))) { \
              results.push_back(std::make_pair(#NAME, show?0: NAME::all_tests(echo))); \
          }
          add_module_test(warnings);
          add_module_test(control);
          add_module_test(data_view);
          add_module_test(icomap);
          add_module_test(icogrid);
          add_module_test(bitmap);
          add_module_test(color);
          add_module_test(impera);
          add_module_test(window);
#undef    add_module_test
      } // testing scope

      int status(0);
      if (results.size() < 1) { // nothing has been tested
          error("test for '%s' not found, use -t '?' to show available module names\n", module);
          status = -1;
      } else {
          if (echo > 0) printf("\n\n#%3ld modules %s tested:\n", results.size(), show?"can be":"have been");
          for(auto result : results) {
              auto const stat = result.second;
              if (echo > 0) {
                  if (show) printf("#    module= %s\n", result.first);
                  else      printf("#    module= %-24s status= %i\n", result.first, int(stat));
              } // echo
              status += std::abs(int(stat));
          } // result
          if (echo > 0) {
              if (!show) printf("\n#%3ld modules have been tested,  total status= %d\n\n", results.size(), int(status));
              if (status > 0) printf("# Warning! At least one module test failed!\n");
          } // echo
      } // something has been tested
      return status;
  } // run_unit_tests

  int show_help(char const *executable) {
      printf("Usage %s [OPTION]\n"
        "   --help           [-h]\tThis help message\n"
        "   --test <module>  [-t]\tTest module, e.g. impera.\n"
        "   --verbose        [-v]\tIncrement verbosity level\n"
        "   +<name>=<value>      \tModify variable environment\n"
        "\n", executable);
      return 0;
  } // show_help

  int show_version(char const *executable="#", int const echo=0) {
#ifdef _GIT_KEY
      // stringify the value of a macro, two expansion levels needed
      #define macro2string(a) stringify(a)
      #define stringify(b) #b
      auto const git_key = macro2string(_GIT_KEY);
      #undef  stringify
      #undef  macro2string
      control::set("git.key", git_key); // store in the global variable environment
      if (echo > 0) std::printf("# %s git checkout %s\n\n", executable, git_key);
#endif // _GIT_KEY
      return 0;
  } // show_version

} // namespace main_tools

int main(int const argc, char const *argv[]) {
    status_t stat(0);
    char const *test_unit = nullptr; // the name of the unit to be tested
    int run_tests{0};
    int verbosity{3}; // set low
    if (argc < 2) {
        printf("%s: no arguments passed! try --help\n", (argc < 1)?__FILE__:argv[0]); 
        return -1;
    } // no argument passed to executable
    for(int iarg = 1; iarg < argc; ++iarg) {
        assert(nullptr != argv[iarg]);
        char const ci0 = *argv[iarg]; // char #0 of command line argument #1
        if ('-' == ci0) {

            // options (short or long)
            char const ci1 = *(argv[iarg] + 1); // char #1 of command line argument #1
            char const IgnoreCase = 32; // use with | to convert upper case chars into lower case chars
            if ('-' == ci1) {

                // long options
                std::string option(argv[iarg] + 2); // remove two '-' in front
                if ("help" == option) {
                    return main_tools::show_help(argv[0]);
                } else 
                if ("version" == option) {
                    return main_tools::show_version(argv[0], 1);
                } else 
                if ("verbose" == option) {
                    verbosity = 6; // set high
                } else
                if ("test" == option) {
                    ++run_tests; if (iarg + 1 < argc) test_unit = argv[iarg + 1];
                } else {
                    ++stat; warn("# ignored unknown command line option --%s", option.c_str());
                } // option

            } else { // ci1

                // short options
                if ('h' == (ci1 | IgnoreCase)) {
                    return  main_tools::show_help(argv[0]);
                } else
                if ('v' == (ci1 | IgnoreCase)) {
                    ++verbosity; verbosity += 3*('V' == ci1); // increment by 'V':4, 'v':1
                } else
                if ('t' == (ci1 | IgnoreCase)) {
                    ++run_tests; if (iarg + 1 < argc) test_unit = argv[iarg + 1];
                } else {
                    ++stat; warn("# ignored unknown command line option -%c", ci1);
                } // ci1

            } // ci1

        } else // ci0
        if ('+' == ci0) {
            stat += control::command_line_interface(argv[iarg] + 1, iarg); // start after the '+' char
        } else
        if (argv[iarg] != test_unit) {
            ++stat; warn("# ignored command line argument \'%s\'", argv[iarg]);
        } // ci0

    } // iarg
    int const echo = control::get("verbosity", double(verbosity)); // define default verbosity here
    if (echo > 0) {
        printf("\n#");
        for(int iarg = 0; iarg < argc; ++iarg) {
            printf(" %s", argv[iarg]); // repeat the command line arguments
        }   printf("\n");
    } // echo
    if (echo > 0) main_tools::show_version(argv[0], 0); // silent to put the git key into the variable environment
    if (echo > 0) printf("\n# verbosity = %d\n", echo);

    if (run_tests) stat += main_tools::run_unit_tests(test_unit, echo);

    int const control_show = control::get("control.show", 1.); // 0:show none, 1:show used, 2:show unused, 4:show defaults, 6:show all
    if (control_show && echo > 0) control::show_variables(control_show);

    if (echo > 0) warnings::show_warnings(3);
    warnings::clear_warnings(1);

    return int(stat);
} // main
