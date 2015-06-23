#ifndef _ESPotentialMap_h_
#define _ESPotentialMap_h_

#include "Main.h"
#include "Strings.h"
#include "PotentialMap.h"

 
class ESPotentialMap : public PotentialMap {
  public:
    ESPotentialMap(string dx_filename) : PotentialMap(dx_filename) {
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

      double dU1[3], dU2[3];
      dU1[0] = (data[grid[0]-1][grid[1]][grid[2]] - E) / (delta);
      dU1[1] = (data[grid[0]][grid[1]-1][grid[2]] - E) / (delta);
      dU1[2] = (data[grid[0]][grid[1]][grid[2]-1] - E) / (delta);

      dU2[0] = (E - data[grid[0]+1][grid[1]][grid[2]]) / (delta);
      dU2[1] = (E - data[grid[0]][grid[1]+1][grid[2]]) / (delta);
      dU2[2] = (E - data[grid[0]][grid[1]][grid[2]+1]) / (delta);

      *e = E * q;

      double qDiv2 = q / 2;//wrap it all in
      F->x += qDiv2 * (dU1[0] + dU2[0]);
      F->y += qDiv2 * (dU1[1] + dU2[1]);
      F->z += qDiv2 * (dU1[2] + dU2[2]);

      return true;
    }

};

#endif
