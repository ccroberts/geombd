#ifndef _ESPotentialMap_h_
#define _ESPotentialMap_h_

#include "Main.h"
#include "Strings.h"
#include "BinaryPotentialMap.h"

 
class ESPotentialMap : public BinaryPotentialMap {
  public:
    ESPotentialMap(string bpm_filename, string atomtype) : BinaryPotentialMap(bpm_filename, atomtype) {
    }

    bool potential(vertex *R, double q, double *e) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      *e = q * data[grid[0]][grid[1]][grid[2]];

      return true;
    }


    bool force(vertex *R, vertex *F, double q, double *e) {
      int grid[3];
      double E;
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      E = data[grid[0]][grid[1]][grid[2]];

      double dU[3];
      dU[0] = (data[grid[0]-1][grid[1]][grid[2]] - data[grid[0]+1][grid[1]][grid[2]]) / (2. * delta);
      dU[1] = (data[grid[0]][grid[1]-1][grid[2]] - data[grid[0]][grid[1]+1][grid[2]]) / (2. * delta);
      dU[2] = (data[grid[0]][grid[1]][grid[2]-1] - data[grid[0]][grid[1]][grid[2]+1]) / (2. * delta);

      if(e) *e = E * q;

      F->x += q * dU[0];
      F->y += q * dU[1];
      F->z += q * dU[2];

      return true;
    }

};

#endif
