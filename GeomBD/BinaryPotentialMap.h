#ifndef _BinaryPotentialMap_h_
#define _BinaryPotentialMap_h_

#include "Main.h"
#include "Strings.h"

 
class BinaryPotentialMap {
  friend class Model;

  protected:
    bool header = true,
         body = false,
         footer = false;

    int N[3];
    int Nt;
    double origin[3];
    double delta;
    double ***data;

  public:
    string type;

  public:
    BinaryPotentialMap(string bpm_filename, string atomtype) {
      type = atomtype;

      ifstream fd(bpm_filename.c_str(), ios::in | ios::binary);

      fd.read((char*)&origin, sizeof(double) * 3);
      fd.read((char*)&N, sizeof(int) * 3);
      fd.read((char*)&delta, sizeof(double));
      Nt = N[0] * N[1] * N[2];

      data = (double***)calloc(N[0], sizeof(double**));
      for(int nx=0; nx < N[0]; nx++) {
        data[nx] = (double**)calloc(N[1], sizeof(double*));              
        for(int ny=0; ny < N[1]; ny++) {
          data[nx][ny] = (double*)calloc(N[2], sizeof(double));              
        }
      }

      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          for(int nz=0; nz < N[2]; nz++) {
            fd.read((char*)&data[nx][ny][nz], sizeof(double));
            if(data[nx][ny][nz] > 400.) data[nx][ny][nz] = 400.;
          }
        }
      }

    }

    virtual ~BinaryPotentialMap() {
      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          free(data[nx][ny]);
        }
        free(data[nx]);
      }
      free(data);
    }

    bool coordinateToGrid(vertex *R, int *G) {
      G[0] = (int)floor(((R->x-origin[0])/delta) + 0.5);
      if(G[0] <= 0 or G[0] >= N[0]-1) {
        return false;
      }
      G[1] = (int)floor(((R->y-origin[1])/delta) + 0.5);
      if(G[1] <= 0 or G[1] >= N[1]-1) {
        return false;
      }
      G[2] = (int)floor(((R->z-origin[2])/delta) + 0.5);
      if(G[2] <= 0 or G[2] >= N[2]-1) {
        return false;
      }

      return true;
    }

    bool potential(vertex *R, double *e) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      *e = data[grid[0]][grid[1]][grid[2]];

      return true;
    }

    bool force(vertex *R, vertex *F, double *e) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return false;

      double E = data[grid[0]][grid[1]][grid[2]];

      double dU[3];
      dU[0] = (data[grid[0]-1][grid[1]][grid[2]] - data[grid[0]+1][grid[1]][grid[2]]) / (2. * delta);
      dU[1] = (data[grid[0]][grid[1]-1][grid[2]] - data[grid[0]][grid[1]+1][grid[2]]) / (2. * delta);
      dU[2] = (data[grid[0]][grid[1]][grid[2]-1] - data[grid[0]][grid[1]][grid[2]+1]) / (2. * delta);

      *e = E;

      F->x += dU[0];
      F->y += dU[1];
      F->z += dU[2];

      return true;
    }

    void translate(double dx, double dy, double dz) {
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
    }

};

#endif
