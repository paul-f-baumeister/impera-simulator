#pragma once

  // constants only

  double constexpr year=365.25, per_year=1./year;
  
// double constexpr deathrate = 6.103515625e-05; // each human dies within ~ 2^14 days (44.85 a)
  double constexpr deathrate = 4e-5; // each human dies within ~ 3*2^13 days (67.3 a)
  double constexpr birthrate = deathrate * 2;
  double constexpr diffusion = 3.3e-4;
  double constexpr powerlimit = 1e3;
  double constexpr brutality = 1e-6; // must be 10^(-4) or smaller
