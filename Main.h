#ifndef Main_h_
#define Main_h_

#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> 
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <mkl.h>
#include <limits.h>

using namespace std;


/*
 * Constants
 */
#define kB 1.9858775e-3 // kcal/molK
#define kC 332.0623183 // coulomb constant
#define Na 6.0221415e23
#define LperA3 1.0e-27
#define nMperM 1.0e9


/*
 * Conversions
 */
#define CONCENTRATION(Nparticles, BOXL) (((()Nparticles / Na) / ((BOXL * BOXL * BOXL) * LperA3)) * nMperM)
#define BOXSIZEFORnMCONC(nMCONC) pow((nMperM * (1.0 / Na) / (nMCONC * LperA3)), 1./3.)

/*
 * Convenience
 */

inline double random(double rangeStart, double rangeEnd) {
  return (((double)rand() / (double)INT_MAX) * (rangeEnd - rangeStart)) + rangeStart;
}


/*
 * Type definitions
 */

struct vertex {
  double x;
  double y;
  double z;
};

enum SimulationConfig {
  CONFIGURATION_RADIAL,
  CONFIGURATION_ABSOLUTE_RADIAL,
  CONFIGURATION_ABSOLUTE_PERIODIC
};

struct BindingSite {
  double r2;
  double x, y, z;
};



typedef unsigned int unt;

#endif
