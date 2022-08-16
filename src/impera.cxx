#include <cstdio> // std::printf, ::snprintf
#include <cassert> // assert
#include <cstdint> // int64_t, uint32_t, uint8_t
#include <fstream> // std::fstream
#include <sstream> // std::sstream
#include <string> // std::string
#include <vector> // std::vector<T>

#ifdef _OPENMP
    #include <omp.h> // omp_get_num_threads, omp_get_max_threads, omp_get_thread_num
#endif // _OPENMP
#include "mpi_replacements.hxx" // MPI_Wtime, ...

#include "icomap.hxx" // ::create_world_map, GridCell_t, ::n_ico_vertices
#include "config.hxx" // constants ...
#include "bitmap.hxx" // ::write_bmp_file
#include "color.hxx" // ::def, ::colorchar, ::display
#include "data_view.hxx" // view2D<T>, view3D<T>

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "warnings.hxx" // warn
#include "political.hxx" // ::init_political_colors
#include "window.hxx" // ::display3D, ::init, ::finalize

namespace impera {

  using namespace mpi_replacements; // MPI_*

  void display_in_terminal(view3D<float> const & pop_map, int const height, int const width) {
      std::fflush(stdout);
      for (int ih = 0; ih < height; ++ih) {
          std::printf("~~~~");
          for (int iw = 0; iw < width; ++iw) {
              float const red   = pop_map(ih, iw, 2);
              float const green = pop_map(ih, iw, 1);
              float const blue  = pop_map(ih, iw, 0);
              std::printf("%s ", color::colorchar(red, green, blue).c_str());
          } // iw
          std::printf("%s~~~~\n", color::def);
      } // ih
      std::fflush(stdout);
  } // display_in_terminal

  status_t read_resource_file(std::vector<GridCell_t> & resources, char const *filename, int const echo) {
      std::ifstream infile(filename, std::ifstream::in);
      if (infile.fail()) {
          warn("Unable to open file '%s' for reading resources", filename);
          return 1;
      } // failed

      size_t nr{0}; char comment;
      std::string line;
      if(!std::getline(infile, line)) { 
          warn("Unable to read 1st line in file '%s'", filename);
          return 1;
      }
      std::istringstream iss(line);
      if (!(iss >> comment >> nr)) {
          warn("Unable to read comment line in file '%s'", filename);
          return 2;
      }
      if (echo > 0) std::printf("# try to read %ld regions from file '%s'\n", nr, filename); 
      size_t const nregions = nr;
      resources.resize(nregions);
      size_t linenumber{0};
      int constexpr nn_ALL = 7;
      size_t nnhist[1 + nn_ALL] = {0,0,0,0, 0,0,0,0};
      nr = 0; // reset region counter
      while (std::getline(infile, line) && (nr < nregions)) {
          ++linenumber;
          int64_t gid, n[6]; double r;
          std::istringstream iss(line);
          if (!(iss >> gid >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5] >> r)) {
              std::printf("# failed parsing in %s:%d reads \"%s\", stop\n", filename, linenumber, line.c_str());
              break; // error
          } else {
              if (gid > -3) {
                  if (gid < 0) gid += nregions; // NorthPole and SouthPole
                  if (gid < nregions) {
                      assert(gid >= 0);
                      assert(gid < resources.size());
                      auto & cell = resources[gid];
                      cell.resources = std::max(0.0, r);
                      int nn{0};
                      for (int i6 = 0; i6 < 6; ++i6) {
                          cell.neighbors[i6] = 0; // init
                          if (n[i6] > -3) {
                              cell.neighbors[nn] = (n[i6] + nregions) % nregions;
                              ++nn;
                          } // valid index
                      } // i6
                      cell.n_neighbors = nn;
                      ++nnhist[nn]; ++nnhist[nn_ALL];

                      ++nr; // a cell has been successfully initialized
                  } else {
                      warn("GridCell index in %s:%d out of bounds [-2, %ld) found %ld", filename, linenumber, nregions, gid);
                  }
              } else {
                  warn("GridCell index in %s:%d out of bounds [-2, %ld) found %ld", filename, linenumber, nregions, gid);
              } // gid > -3
          } // line parsing with iss failed
      } // while

      // show nnhist
      if (echo > 2) {
          std::printf("# found "); 
          size_t nnall{0};
          for (int nn = 0; nn < nn_ALL; ++nn) {
              nnall += nnhist[nn];
              if (nnhist[nn] > 0) std::printf("%d neighbors (%ld times), ", nn, nnhist[nn]);
          } // nn
          std::printf("total=%ld in '%s'\n", nnall, filename); 
          assert(nnall == nnhist[nn_ALL]);
      } // echo

      return nr - nregions; // success if the expected number of regions has been found
  } // read_resource_file




  template <typename real_t, int Nspecies=3>
  class Impera_t {
  private:
      int RenderTime, BitMapTime, PoliticalMapTime, DisplayTime, DisplayPopTime, DisplayPerformance;
      char const *BitMapFileName, *PicturePath, *ResourceFileBaseName, *ResourceFileExtension, *ResourceFilePath;
      bool DisplayPolitical;
      unsigned Level; // icosahedral grid level
      static real_t constexpr tax = 0.002; // equilibration of power (per day tax)  --> if we assume that we have to pass 33% per year, this is 0.66**(1./365.25) = 1-0.001137
      std::vector<GridCell_t> resources;
      size_t Nregions;
      // Quantities with (Nspecies,Nregions)
      view2D<real_t> pop;  // population(s,r)
      view2D<real_t> pop1; // population(s,r)
      view2D<real_t> pot;  // potential(s,r)
      view2D<real_t> pot1; // potential(s,r)
      view2D<real_t> pwr;  // power(s,r)
      view2D<real_t> mobil; // mobilization(s,r)
      // Quantities with (Nspecies,Nspecies)
      view2D<double> grad_pop;  // (s,s)
      view2D<double> grad_pwr;  // (s,s)
      view2D<double> hostility; // (s,s)
      view2D<double> hostile;   // (s,s)
      view2D<double> casulties; // (s,s)
      view2D<float>  political_colors; // (Nspecies + 1, 4);
      // Quantities with (Nspecies)
      std::vector<double> power_at_spawn; // (s)
      std::vector<double> allpower; // (s)
      std::vector<double> allpopul; // (s)
      std::vector<double> perpopul; // (s)

  public:

      Impera_t(int const echo=1, int const run=1) // constructor
          : RenderTime           (control::get("RenderTime", -1.)) // 32:every month, 0:never, 1: every day, 356: every year
          , BitMapTime           (control::get("BitMapTime", 1461.)) // 365*4+1 == 4 years exactly
          , PoliticalMapTime     (control::get("PoliticalMapTime", 1461.)) // 365*4+1 == 4 years exactly, 0:never
          , DisplayTime          (control::get("DisplayTime", 32.)) // 32:every month, 0:never, 1: every day, 356: every year
          , DisplayPopTime       (control::get("DisplayPopTime", 32.)) // 32:every month, 0:never, 1: every day, 356: every year
          , DisplayPerformance   (control::get("DisplayPerformance", 255.))
          , BitMapFileName       (control::get("BitMapFileName", "img"))
          , PicturePath          (control::get("PicturePath", "../pic"))
          , ResourceFileBaseName (control::get("ResourceFileBaseName", "mapfile_L"))
          , ResourceFileExtension(control::get("ResourceFileExtension", ".dat"))
          , ResourceFilePath     (control::get("ResourceFilePath", "../icogrid"))
          , DisplayPolitical     (control::get("DisplayPolitical", 0.)) // 0:never, 1:always
          , Level                (control::get("Level", 3.)) // icosahedral grid level
    { // constructor body

        status_t stat(0);
        char infilename[96];
        std::snprintf(infilename, 95, "%s/%s%i%s", ResourceFilePath, ResourceFileBaseName, Level, ResourceFileExtension);
        if (echo > 0) std::printf("# read resource map from '%s'\n", infilename);

        stat += read_resource_file(resources, infilename, echo);
        size_t const nregions = resources.size();

        { // scope: VerifyResourceMap
            int const VerifyResourceMap = control::get("VerifyResourceMap", 1.);
            if (echo > 0) std::printf("# VerifyResourceMap=%d\n", VerifyResourceMap);
            if (VerifyResourceMap) {
                assert(nregions == icomap::n_ico_vertices(Level));
                int const h = icomap::map_height(Level),
                          w = icomap::map_width(Level);
                view2D<float> resource_ico(nregions, 4, 0.f);
                double scale{0};
                bool const output_resource = (control::get("OutputResource", 0.) > 0);
                for (size_t ir = 0; ir < nregions; ++ir) {
                    double const r = std::sqrt(std::sqrt(resources[ir].resources));
                    resource_ico(ir,0) = r; // red
                    resource_ico(ir,1) = 0; // green
                    resource_ico(ir,2) = 0; // blue
                    if (output_resource) std::printf("%g\n", r);
                    scale = std::max(scale, r);
                } // ir
                std::printf("# largest resource value is %g\n", pow4(scale));

                int const nhist = control::get("Histogram", 32.);
                if (nhist > 1) {
                    std::vector<int> hist(nhist, 0);
                    double const factor = (nhist - .99)/scale;
                    for (size_t ir = 0; ir < nregions; ++ir) {
                        double const r = std::sqrt(std::sqrt(resources[ir].resources));
                        int const ihist = r*factor;
                        ++hist[ihist];
                    } // ir
                    double const percent = 100./nregions, invfactor = scale/(nhist - 1.);
                    std::printf("# histogram of resources (%d bins) in %%:\n", nhist);
                    for (int ihist = 0; ihist < nhist; ++ihist) {
                        std::printf("%g %g\n", pow4(ihist*invfactor), hist[ihist]*percent);
                    } // ihist
                } // create a histogram

                float const inv_scale = 1./scale;
                view3D<float> resource_map(h, w, 4, 0.f);
                bool const icosahedral = true;
                stat += icomap::create_world_map(resource_map, resource_ico, Level, 1.f, inv_scale, icosahedral);
                char outfilename[96];
                std::snprintf(outfilename, 95, "%s/ico_resource-L%i", PicturePath, Level);
                stat += bitmap::write_bmp_file(outfilename, resource_map.data(), h, w);

                stat += icomap::create_world_map(resource_map, resource_ico, Level, 1.f, inv_scale);
                std::snprintf(outfilename, 95, "%s/ico_resources_square-L%i", PicturePath, Level);
                stat += bitmap::write_bmp_file(outfilename, resource_map.data(), h, w);

                display_in_terminal(resource_map, h, w);

                if (VerifyResourceMap < 0) return;
            } // VerifyResourceMap
        } // scope

#ifdef _OPENMP
        #pragma omp master
            std::printf("# running with max %d threads\n", omp_get_max_threads());
#endif // _OPENMP

#ifdef _MPI
          MPI_Comm const comm = MPI_COMM_WORLD;
          int MPI_communicator_size{1}, MPI_communicator_rank{0};
          MPI_Comm_size(comm, &MPI_communicator_size);
          MPI_Comm_rank(comm, &MPI_communicator_rank);

          int const MPIsize = std::max(1, MPI_communicator_size);
          int const MPIrank = MPI_communicator_rank;
#else
          int constexpr MPIsize = 1;
#endif

          // MPI distribute regions
          Nregions = (nregions - 1)/MPIsize + 1;
          size_t const Nregions_aligned = align<5>(Nregions); // multiples of 32

          // main memory consumers
          pop   = view2D<real_t>(Nspecies,Nregions_aligned, real_t(1e-6)); // population(s,r), init as very sparsely populated
          pop1  = view2D<real_t>(Nspecies,Nregions_aligned); // population(s,r)
          pot   = view2D<real_t>(Nspecies,Nregions_aligned, real_t(0)); // potential(s,r)
          pot1  = view2D<real_t>(Nspecies,Nregions_aligned); // potential(s,r)
          pwr   = view2D<real_t>(Nspecies,Nregions_aligned, real_t(1)); // power(s,r)
          mobil = view2D<real_t>(Nspecies,Nregions_aligned);  // mobilization(s,r)

          int const Nspecies_aligned = align<0>(Nspecies); // ToDo: use Nspecies_aligned in the lower allocations
          grad_pop  = view2D<double>(Nspecies,Nspecies_aligned);
          grad_pwr  = view2D<double>(Nspecies,Nspecies_aligned);
          hostility = view2D<double>(Nspecies,Nspecies_aligned, 0.0);
          hostile   = view2D<double>(Nspecies,Nspecies_aligned);
          casulties = view2D<double>(Nspecies,Nspecies_aligned, 0.0); // (s,s)

          // std::vector<double> share(Nspecies) ; //, allpower, allpopul, perpopul, strength, weights, relative_force, tax, power_at_spawn=0. ! (s)
          power_at_spawn = std::vector<double>(Nspecies, 0.0); // (s)
          allpower       = std::vector<double>(Nspecies, 0.0); // (s)
          allpopul       = std::vector<double>(Nspecies, 0.0); // (s)
          perpopul       = std::vector<double>(Nspecies, 0.0); // (s)


          // colors for the political map
          political_colors = view2D<float>(Nspecies + 1, 4, 0.f);
          assert(0 < political::init_political_colors(political_colors.data(), Nspecies));

          double const RescaleResources = control::get("RescaleResources", 1.0);
          if (1 != RescaleResources) {
              double sum{0};
              #pragma omp parallel for reduction(+:sum)
              for(size_t ir = 0; ir < Nregions; ++ir) {
                  resources[ir].resources *= RescaleResources;
                  sum += resources[ir].resources;
              } // ir
              std::printf("# sum of all resources (rescaled) %g\n", sum);
          } // RescaleResources


          { // scope: initialize
              float start_pop{1e3}; float t{9e37};
              float const divisor = (Nspecies > 12) ? 0.75f : 0.5f; // fails to place 32 species with 0.5
              for (int is = 0; is < Nspecies; ++is) { // serial loop over species, loop-carried dependency
                  int ivr{-1}; float maxres{0};
                  #pragma omp parallel for
                  for(size_t ir = 0; ir < Nregions; ++ir) {
                      auto const res = resources[ir].resources;
                      #pragma omp critical
                      if (res < t && res > maxres) { maxres = res; ivr = ir; }
                      pwr(is,ir) = is ? pwr(is - 1,ir)/2 : 100;
                  } // ir
                  assert(ivr > -1);
                  pop(is,ivr) = start_pop; // start in a fruitful region

                  std::printf("# population #%i starts in region #%i with %g resources\n", is, ivr, maxres);

                  t = resources[ivr].resources * divisor;
                  start_pop *= divisor;
              } // is
          } // scope

//           if (RenderTime > 0) {
//               std::printf("# window::init()\n");
//               window::init();
//           } // render_window

          if (run) {
              run_simulation();
          } // run

      } // constructor

      status_t run_simulation() {
          status_t stat(0);
          double const nThousandYears = control::get("nThousandYears", 1.);
          uint64_t const number_of_days = year*1000*nThousandYears;

          double time_all_start = MPI_Wtime();

          for(uint64_t it = 0; it < number_of_days; ++it) {
              stat += time_loop(it, &time_all_start);
          } // it

          auto const average_speed = number_of_days/(MPI_Wtime() - time_all_start);
          #pragma omp master
          std::printf("\n# Average performance %g simulated days per second\n", average_speed);

          return finalize();
      } // run_simulation


      status_t finalize() {
          status_t stat(0);
          // show hostility and casulties
          for (int hc = 0; hc < 2; ++hc) {
              std::printf("\n# %s\n", hc?"casulties":"hostility");
              for (int is = 0; is < Nspecies; ++is) {
                  std::printf("%6i", is);
                  double csum{0};
                  for (int js = 0; js < Nspecies; ++js) {
                      if (hc) std::printf("%12.3e", casulties(js,is));
                      else    std::printf("%12.6f", hostility(js,is));
                      csum += casulties(js,is);
                  } // js
                  if (hc) std::printf("\t%16.3e", csum); // sum of casulties
                  std::printf("\n");
              } // is
          } // hc

//           if (RenderTime > 0) {
//               std::printf("# window::finalize()\n");
//               window::finalize();
//           } // render_window

          std::printf("\n# Level = %d done\n", Level);
          return stat;
      } // finalize




  status_t time_loop(
        uint64_t const it // day counter
      , double const *time_all_start=nullptr
  ) { // time-loop, count in days 
      status_t stat(0);

      auto const time_iter_start = MPI_Wtime();

      set(hostile, Nspecies, 0.0); // init day change

      { // begin omp parallel (in the future)


      #pragma omp parallel for
      for (size_t ir = 0; ir < Nregions; ++ir) {

          std::vector<double> strength(Nspecies), weights(Nspecies), relative_force(Nspecies), share(Nspecies);
          double sumstrength{0}, local_pop{0};
          for (int is = 0; is < Nspecies; ++is) {
              strength[is] = pop(is,ir)*pwr(is,ir);
              sumstrength += strength[is];
              local_pop += pop(is,ir);
          } // is
          double const t = (sumstrength > 0) ? 1./sumstrength : 0;
          for (int is = 0; is < Nspecies; ++is) {
              weights[is] = strength[is] * t;
              relative_force[is] = pwr(is,ir) * t;
          } // is

          double local_harvest{0}, d_harvest_d_population{0};
          double const res = resources[ir].resources;
          if (res > 0) {
              local_harvest          = 2.0 *       res / (res + local_pop) * local_pop; // 2RP/(R+P)
              d_harvest_d_population = 2.0 * pow2( res / (res + local_pop) );
          } // resources > 0

          set(share.data(), Nspecies, weights.data(), local_harvest);

          for (int is = 0; is < Nspecies; ++is) {
              mobil(is,ir) = 0;
          } // is

          for (int is = 0; is < Nspecies; ++is) {
              for (int js = 0; js < Nspecies; ++js) {
                  grad_pop(js,is) = weights[is] * ( d_harvest_d_population - relative_force[js] * local_harvest );
                  grad_pwr(js,is) = - local_harvest * ( weights[is] * pop(js,ir) * t );
              } // js
              grad_pwr(is,is) += local_harvest * pop(is,ir) * t;
              grad_pop(is,is) += local_harvest * pow2(relative_force[is]);
              for (int js = 0; js < Nspecies; ++js) {
                  hostile(js,is) -= grad_pop(js,is) * pop(is,ir); // reduction over regions weighted with relative population
                  mobil(js,ir) += grad_pop(js,is) * strength[is]; // weighted with local absolute strength
              } // js
          } // is

          for (int is = 0; is < Nspecies; ++is) {
              auto const surv = std::min(double(pop(is,ir)), share[is]); // survive
              share[is] -= surv; // what is left of the share can be used to increase the power
              auto const birth = birthrate * (1. - pow2((pwr(is,ir) - power_at_spawn[is])/powerlimit));
              auto const per_capita = (pop(is,ir) > 0.)? 1./pop(is,ir) : 0;
              auto const hunger = pop(is,ir) - surv;
              auto const death = std::min(1., deathrate + hunger * per_capita);

              pwr(is,ir) += std::max(0., share[is] * per_capita * 0.001); // change the power
              pop(is,ir) = std::max(0., pop(is,ir) * ((1. - death) * (1. + birth))); // change population
          } // is

      } // ir

    // pot1 = 0.75 * pot + 0.25 * mobil // slowly adopt new gradients
//  #pragma omp parallel
    {
        set(pot1.data(), Nspecies*pot1.stride(), pot.data(), real_t(0.75));
        add_product(pot1.data(), Nspecies*pot1.stride(), mobil.data(), real_t(0.25));
    }

    set(pot, Nspecies, real_t(0));
    real_t const cn_nn[8] = {1., .5, 1/3., .25, .2, 1/6., 1/7.};
    for (size_t ir = 0; ir < Nregions; ++ir) {
        int const nn = resources[ir].n_neighbors;
//      real_t const cn = 1./(nn + 2.), c0 = 2*cn; // normalized: c0 + nn*cn == 1.0
        real_t const c0 = 0.25, cn = 0.75*cn_nn[nn]; // normalized: c0 + nn*cn == 1.0

        for (int is = 0; is < Nspecies; ++is) {
            pot(is,ir) += c0 * pot1(is,ir); // diffusion central
        } // is
        for (int in = 0; in < nn; ++in) {
            auto const jr = resources[ir].neighbors[in];
            for (int is = 0; is < Nspecies; ++is) {
                pot(is,jr) += cn * pot1(is,ir); // diffusion hopping
            } // is
        } // in
    } // ir

    set(pop1.data(), Nspecies*pop1.stride(), pop.data(), real_t(1. - diffusion)); // those who stay

    for (size_t ir = 0; ir < Nregions; ++ir) {
        int const nn = resources[ir].n_neighbors;
        assert(nn <= 6);
        for (int is = 0; is < Nspecies; ++is) {
            double diff[6]; // ToDo: 6:max_N_neighbors
            double sumdiff{0};
            for (int in = 0; in < nn; ++in) {
                auto const jr = resources[ir].neighbors[in]; // target region
                auto const grad = pot(is,jr) - pot(is,ir); // gradient of the potential
                diff[in] = 0.5 + 0.5*grad/(1. + std::abs(grad)); // cheap logistic function
                sumdiff += diff[in];
            } // in
            auto const t = diffusion/sumdiff; // normalize for each species
            for (int in = 0; in < nn; ++in) {
                diff[in] *= t; // scale
                auto const jr = resources[ir].neighbors[in]; // target region
                pop1(is,jr) += pop(is,ir) * diff[in];
            } // in
        } // is
    } // ir

    std::swap(pop, pop1); // adopt diffusion

    // init global sums
    for (int is = 0; is < Nspecies; ++is) {
        allpopul[is] = 0;
        allpower[is] = 0;
    } // is
    double maxpop{0};

    //============WAR===============================================================
    for (size_t ir = 0; ir < Nregions; ++ir) {
        for (int is = 0; is < Nspecies; ++is) { // loop over aggressors
            auto const si = pop(is,ir)*pwr(is,ir); // strength of the aggressor species is
            if (si >= 1e-9) {
                auto const t = 1./si;
                for (int js = 0; js < Nspecies; ++js) {  // loop over victims
                    auto const war_mercy = (js == is) ? 0.125 : 1.0;
                    auto const sj = pop(js,ir)*pwr(js,ir); // strength of species js
                    if (sj >= 1e-9) {
                        auto const war_efficiency = 1./(1. + pow2(std::log( sj*t ))); // this is symmetric under exchange is <--> js
                        auto const war_intensity = std::max(0.0, hostility(js,is)); // the effect of polarizable media ...
                        auto const x = pwr(is,ir)/pwr(js,ir);
                        auto const saturation = x/(std::abs(x) + 1.);
                        double war_casulties = brutality * saturation * pop(js,ir) * pop(is,ir) * war_efficiency * war_intensity * war_mercy;
                        war_casulties = std::min(war_casulties, double(pop(js,ir)));

                        pop(js,ir) = std::max(0., pop(js,ir) - war_casulties); // people get killed
                        casulties(is,js) += war_casulties;
                    } // sj
                } // js
            } // si
        } // is

        for (int is = 0; is < Nspecies; ++is) {        
            allpopul[is] += pop(is,ir);
            allpower[is] += pop(is,ir) * pwr(is,ir);
            maxpop = std::max(maxpop, double(pop(is,ir)));
        } // is
    } // ir
    // ============WAR===============================================================


    for (int is = 0; is < Nspecies; ++is) {       
        perpopul[is] = (allpopul[is] > 0) ? 1./allpopul[is] : 0;
        allpower[is] *= perpopul[is]; // normalized
        for (int js = 0; js < Nspecies; ++js) {
            hostile(js,is) *= perpopul[is]; // ! average over the population
        } // js

        for (size_t ir = 0; ir < Nregions; ++ir) {
            // step-wise equilibration of power
            pwr(is,ir) = (1 - tax) * pwr(is,ir) + tax * allpower[is];
        } // ir
    } // is



    bool const display_screen = (DisplayTime > 0) && (0 == (it % DisplayTime));
    bool const display_pop = (DisplayPopTime > 0) && (0 == (it % DisplayPopTime));
    if (display_screen || display_pop) {
        std::printf("# year %.2f pop ", it*per_year);
        double world_pop{0};
        for (int is = 0; is < Nspecies; ++is) {        
            std::printf(" %.2e", allpopul[is]);
            bool constexpr show_power_with_colors = true;
            if (show_power_with_colors) {
                auto const ccc = color::colorchar(0 == (is%3), 0 == ((is + 1)%3), 0 == ((is + 2)%3));
                std::printf("%s@%.1f%s", ccc.c_str(), allpower[is], color::def);
            } // show_power_with_colors
            world_pop += allpopul[is];
        } // is
        std::printf("%s sum= %.4e%s\n", color::def, world_pop, color::def);
    } // display_screen or display_pop

    bool const export_bitmap = (BitMapTime > 0) && (0 == (it % BitMapTime));
    bool const render_window = (RenderTime > 0) && (0 == (it % RenderTime));
    if (export_bitmap || display_screen || render_window) {

        std::vector<uint8_t> color_mapping(Nspecies);
        for (int is = 0; is < Nspecies; ++is) {
            color_mapping[is] = is % 3;
        } // is
        view2D<float> pop_ico(Nregions, 3, 0.f);
        float max_p{0};
        for (size_t ir = 0; ir < Nregions; ++ir) {
            auto const pi = pop_ico[ir];
            for (int is = 0; is < Nspecies; ++is) {
                float const p = std::sqrt(std::sqrt(pop(is,ir)));
                pi[color_mapping[is]] += p;
            } // is
            max_p = std::max(std::max(std::max(pi[0], pi[1]), pi[2]), max_p);
        } // is
        float const scale_p = 1./max_p;
        int const h = icomap::map_height(Level),
                  w = icomap::map_width(Level);
        view3D<float> pop_map(h, w, 4, 0.f);
        stat += icomap::create_world_map(pop_map, pop_ico, Level, 1.f, scale_p);

        if (export_bitmap) {

            char bitmapfile[512];
            std::snprintf(bitmapfile, 512, "%s/%s-L%i.%06li", PicturePath, BitMapFileName, Level, it/BitMapTime);
            stat += bitmap::write_bmp_file(bitmapfile, pop_map.data(), h, w);

        } // export bitmap

        if (display_screen) {

            display_in_terminal(pop_map, h, w);

        } // display_screen

        if (render_window) {
            if (icomap::n_ico_vertices(Level) == Nregions) {
//              std::printf("# window::display3D(Level, pop_ico.data());\n"); // DEBUG
                window::display3D(Level, pop_ico.data());
            } else warn("Rendering only supported for icosahedral grids");
        } // render in OpenGL window

    } // export bitmap or show in terminal or render

    bool const export_political = ((PoliticalMapTime > 0) && (0 == (it % PoliticalMapTime)));
    if (export_political) {
        // construct an overlap matrix measuring if two species are next to each other or living on the same lands
        view2D<double> ovl(Nspecies, Nspecies, 0.0);
        for (size_t ir = 0; ir < Nregions; ++ir) {
            for (int is = 0; is < Nspecies; ++is) {
                double const pop_is_ir = pop(is,ir);
                for (int js = 0; js < Nspecies; ++js) {
                    ovl(is,js) += pop_is_ir*pop(js,ir);
                } // js
            } // is
        } // ir

        std::vector<double> factor(Nspecies, 1.0);
        for (int is = 0; is < Nspecies; ++is) {
            factor[is] = 1./std::sqrt(ovl(is,is));
        } // is

        // renormalize to have unity on the diagonal
        for (int is = 0; is < Nspecies; ++is) {
            for (int js = 0; js < Nspecies; ++js) {
                ovl(is,js) *= factor[is]*factor[js];
            } // js
        } // is

        if (1) {
            std::printf("# year %.2f overlap (*10):\n", it*per_year);
            for (int is = 1; is < Nspecies; ++is) {
                std::printf("#%3i  ", is);
                for (int js = 0; js < is; ++js) {
                    std::printf("%6.3f", 10*ovl(is,js)); // show only off-diagonal elements, matrix is symmetric and has 1s on the diagonal
                } // js
                std::printf("\n");
            } // is
        } // output to screen

        assert(0 < political::update_political_colors(political_colors.data(), Nspecies, ovl.data()));

        // find out which cell has which color in the political map
        std::vector<int16_t> strongest_species(Nregions, -1); // -1:uninitialize
        for (size_t ir = 0; ir < Nregions; ++ir) {
            real_t max_strength{1e-9}; // increase this threshold to slow down the conquest of uninhabitated land (and sea).
            int strongest{Nspecies}; // Nspecies:none
            if (resources[ir].resources > 0) {
                for (int is = 0; is < Nspecies; ++is) {
                    real_t const strength = pop(is,ir)*pwr(is,ir);
                    if (strength > max_strength) {
                        max_strength = strength;
                        strongest = is;
                    } // stronger
                } // is
            } // resources > 0
            strongest_species[ir] = strongest;
        } // ir

        // now pick Nspecies + 1 different colors to color the map (last color for unexplored terrain);
        // these RGB vector can move in the color cube but forbidden are white and black by repulsion
        // and also there is a stronger repulsion of two species that have an overlap.
        view2D<float> pop_ico(Nregions, 4, 0.f);
        for (size_t ir = 0; ir < Nregions; ++ir) {
            set(pop_ico[ir], 3, political_colors[strongest_species[ir]]);
        } // ir
        int const h = icomap::map_height(Level), w = icomap::map_width(Level);
        view3D<float> pol_map(h, w, 4, 0.f);
        stat += icomap::create_world_map(pol_map, pop_ico, Level, 1.f);

        char bitmapfile[512];
        std::snprintf(bitmapfile, 512, "%s/political-L%i.%06li", PicturePath, Level, it/PoliticalMapTime);
        stat += bitmap::write_bmp_file(bitmapfile, pol_map.data(), h, w);

        if (DisplayPolitical > 0) {
            display_in_terminal(pol_map, h, w);
        } // DisplayPolitical

    } // export_political

    scale(hostility.data(), Nspecies*Nspecies, 0.9); // here, the hysteresis comes in.
    add_product(hostility.data(), Nspecies*Nspecies, hostile.data(), 0.1);


//!+respawninlined
    {
      real_t constexpr population_threshold = 1;

      int ks{-1}; // impossible value (-1:not found)
      for (int is = 0; is < Nspecies; ++is) {
          if (allpopul[is] < population_threshold) ks = is;
      } // is
      if (ks > -1) {

          int jr{-1}; // impossible value (-1:not found)
          float maxval_local_pop{0};
          for (size_t ir = 0; ir < Nregions; ++ir) {
              double local_pop{0};
              for (int is = 0; is < Nspecies; ++is) {
                  local_pop += pop(is,ir);
              } // is
              if (local_pop >= maxval_local_pop) {
                  jr = ir; 
                  maxval_local_pop = local_pop;
              }
          } // ir
          if (jr < 0) error("fatal! most or least populated region not found, choose region #%i", jr); 
          // different from impera.F90 which used cpu_time() as random number generator

          // find the most powerful population in this region (which is not ks)
          int js{-1}; // impossible value (-1: not found)
          // respawn as split of the strongest
          float maxval_local_strength{0};
          for (int is = 0; is < Nspecies; ++is) {
              if (is != ks && pwr(is,jr)*pop(is,jr) >= maxval_local_strength) { js = is; maxval_local_strength = pwr(is,jr)*pop(is,jr); }
          } // is
          if (js < 0) error("fatal! most powerful or strongest population not found");

          auto const strength_transfer = pop(js,jr) * pwr(js,jr) * 0.55f;
          auto const population_transfer = pop(js,jr) * 0.45f;

          std::printf("#\n# species #%i is extinct\n# replace by a split of species #%i\n#\n", ks, js);

          pop(ks,jr) = population_transfer; // split their population (pwr(ks,jr) << 1 anyway)
          pwr(ks,jr) = strength_transfer/population_transfer; // copy their power and give a little more power
          pwr(js,jr) = pwr(js,jr) - strength_transfer/pop(js,jr); // take away the power locally (will be equilibrated soon)
          pop(js,jr) = pop(js,jr) - population_transfer; // take away the population
          power_at_spawn[ks] = pwr(ks,jr); // define power offset, so a new species will again be fertile

          // reset also hostilities
          for (int is = 0; is < Nspecies; ++is) {
              hostility(is,ks) = hostility(is,js);   // copy old hostilities from species (j) (take over their memory)
          } // is
          hostility(ks,js) = hostility(js,js); // immediately after splitting, species (j) will regard the new species (k) as if they were their own people   
          hostility(ks,ks) = 0;                // the new species (k) stands together, no self-hostility, which is a bonus
          hostility(js,ks) = hostility(js,js); // the new species (k) hates the old one (j) as much as usual

      } // create a new population

    }
//!-respawninlined

    } // end omp parallel

    auto const time_iter_end = MPI_Wtime();
    if (time_all_start && 0 == (it & DisplayPerformance)) {
        #pragma omp master
        {
            auto const itr_took =  time_iter_end - time_iter_start;
            auto const avg_took = (time_iter_end - *time_all_start)/(it + 1.);
            std::printf("# Iteration %i took %g seconds, average %g seconds, average speed %g days/sec\n", it, itr_took, avg_took, 1./avg_took);
        } // master
    }

    return stat;
  } // it time-loop


  }; // class Impera_t



  template <int NSpecies>
  void* _run1day(int const echo) {
      static Impera_t<double,NSpecies> *simulation = nullptr;
      static uint64_t iday = 0;
      if (nullptr == simulation) { // on first call
          simulation = new Impera_t<double,NSpecies>(echo, 0);
          iday = 0;
          return (void*)simulation;
      } else { // first call
          simulation->time_loop(iday);
          ++iday;
          return nullptr;
      }
  } // run_one_day

  void* run_one_day(int const Nspecies, int const echo) {
      switch (Nspecies) {
          case  9:  return _run1day<9>(echo);
          case  6:  return _run1day<6>(echo);
          case  3:  return _run1day<3>(echo);
          default:  warn("Requested Nspecies= %d but only {3, 6, 9} instanciated", Nspecies);
                    return nullptr;
      } // switch Nspecies
  } // run_one_day


  template <typename real_t, int Nspecies=3>
  status_t start(int const echo=1, int const run=1) {
      Impera_t<real_t, Nspecies> simulation(echo, run);
      return 0;
  } // start

  status_t test_start(int const echo) {
      bool const dp = ('d' == (*control::get("real_t", "double") | 32));
      int const Nspecies = control::get("Nspecies", 3.);
      switch (Nspecies) {
//        case 32:  return dp ? start<double,32>(echo) : start<float,32>(echo);
//        case 12:  return dp ? start<double,12>(echo) : start<float,12>(echo);
//        case  9:  return dp ? start<double, 9>(echo) : start<float, 9>(echo);
//        case  6:  return dp ? start<double, 6>(echo) : start<float, 6>(echo);
          default:  return dp ? start<double, 3>(echo) : start<float, 3>(echo);
      } // switch Nspecies
  } // test_start

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_start(echo);
      return stat;
  } // all_tests

} // namespace impera
