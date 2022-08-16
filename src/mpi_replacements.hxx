#pragma once

#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::min

#ifdef _OPENMP
    #include <omp.h> // omp_get_wtime
#else
    #include <chrono> // std::chrono::...
#endif

namespace mpi_replacements {
  
    typedef int64_t MPI_Comm;
    
    MPI_Comm constexpr MPI_COMM_NULL = 0, MPI_COMM_WORLD = -1;
  
    inline int MPI_Init(int *argc, char ***argv) { return 0; }

    inline int MPI_Finalize() { return 0; }

    inline int MPI_Comm_size(MPI_Comm comm, int *size) { *size = 1; return 0; }
    inline int MPI_Comm_rank(MPI_Comm comm, int *rank) { *rank = 0; return 0; }

    inline double MPI_Wtime() {
#ifdef _OPENMP
        return omp_get_wtime();
#else
        static int init = 0;
        static std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
        if (!init) {
            start_time = std::chrono::high_resolution_clock::now(); // start
        }
        auto const time_now = std::chrono::high_resolution_clock::now(); // stop
        return std::chrono::duration_cast<std::chrono::microseconds>(time_now - start_time).count()*1e-6;
#endif
    } // MPI_Wtime

} // namespace mpi_replacements
