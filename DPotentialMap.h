#ifndef _DPotentialMap_h_
#define _DPotentialMap_h_

#include "Main.h"
#include "Strings.h"
#include "ESPotentialMap.h"

 
class DPotentialMap : public ESPotentialMap {
  public:
    DPotentialMap(string bpm_filename, string atomtype) : ESPotentialMap(bpm_filename, atomtype) {
    }

    bool potential(vertex *R, double q, double *e) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      if(e) *e = fabs(q) * data[grid[0]][grid[1]][grid[2]];

      return true;
    }


    bool force(vertex *R, vertex *F, double q, double *e) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      double E = data[grid[0]][grid[1]][grid[2]];

      double dU[3];
      dU[0] = (data[grid[0]-1][grid[1]][grid[2]] - data[grid[0]+1][grid[1]][grid[2]]) / (2. * delta);
      dU[1] = (data[grid[0]][grid[1]-1][grid[2]] - data[grid[0]][grid[1]+1][grid[2]]) / (2. * delta);
      dU[2] = (data[grid[0]][grid[1]][grid[2]-1] - data[grid[0]][grid[1]][grid[2]+1]) / (2. * delta);

      if(e) *e = E * fabs(q);

      F->x += fabs(q) * dU[0];
      F->y += fabs(q) * dU[1];
      F->z += fabs(q) * dU[2];

      return true;
    }

};

#endif
