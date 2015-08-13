#ifndef _DPotentialMap_h_
#define _DPotentialMap_h_

#include "Main.h"
#include "Strings.h"
#include "ESPotentialMap.h"

 
class DPotentialMap : public ESPotentialMap {
  public:
    DPotentialMap(string bpm_filename) : ESPotentialMap(bpm_filename) {
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

      *e = data[grid[0]][grid[1]][grid[2]];

      double dU1[3], dU2[3];
      dU1[0] = (data[grid[0]-1][grid[1]][grid[2]] - *e) / (delta);
      dU1[1] = (data[grid[0]][grid[1]-1][grid[2]] - *e) / (delta);
      dU1[2] = (data[grid[0]][grid[1]][grid[2]-1] - *e) / (delta);

      dU2[0] = (*e - data[grid[0]+1][grid[1]][grid[2]]) / (delta);
      dU2[1] = (*e - data[grid[0]][grid[1]+1][grid[2]]) / (delta);
      dU2[2] = (*e - data[grid[0]][grid[1]][grid[2]+1]) / (delta);

      if(e) *e *= fabs(q);

      double qDiv2 = fabs(q) / 2;//wrap it all in
      F->x += qDiv2 * (dU1[0] + dU2[0]);
      F->y += qDiv2 * (dU1[1] + dU2[1]);
      F->z += qDiv2 * (dU1[2] + dU2[2]);

      return true;
    }

};

#endif
