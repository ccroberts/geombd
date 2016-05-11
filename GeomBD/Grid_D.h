#ifndef _Grid_D_h_
#define _Grid_D_h_

#include "Main.h"
#include "Strings.h"
#include "Grid.h"

 
class Grid_D : public Grid {
  public:
    Grid_D(string bpm_filename, string atomtype) : Grid(bpm_filename, atomtype) {
    }

    bool potential(vertex *R, double q, double *e) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      *e = fabs(q) * data[grid[0]][grid[1]][grid[2]];

      return true;
    }


    bool force(vertex *R, vertex *F, double q, double *e) {
      int grid[3];
      double E;
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      E = data[grid[0]][grid[1]][grid[2]];

      double dU[3];
      dU[0] = fabs(q) * (data[grid[0]-1][grid[1]][grid[2]] - data[grid[0]+1][grid[1]][grid[2]]) / (2. * delta);
      dU[1] = fabs(q) * (data[grid[0]][grid[1]-1][grid[2]] - data[grid[0]][grid[1]+1][grid[2]]) / (2. * delta);
      dU[2] = fabs(q) * (data[grid[0]][grid[1]][grid[2]-1] - data[grid[0]][grid[1]][grid[2]+1]) / (2. * delta);

      F->x += dU[0];
      F->y += dU[1];
      F->z += dU[2];

      if(e) *e += E * fabs(q);

      return true;
    }

};

#endif
