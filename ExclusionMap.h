#ifndef _ExclusionMap_h_
#define _ExclusionMap_h_

#include "Main.h"
#include "Strings.h"

 
class ExclusionMap {
  protected:
    bool header = true,
         body = false,
         footer = false;

    int N[3];
    int Nt;
    double origin[3];
    double delta;
    short ***data;

  public:
    string type;

  public:
    ExclusionMap(string bpm_filename, string atomtype) {
      type = atomtype;

      ifstream fd(bpm_filename.c_str(), ios::in | ios::binary);

      fd.read((char*)&origin, sizeof(double) * 3);
      fd.read((char*)&N, sizeof(int) * 3);
      fd.read((char*)&delta, sizeof(double));
      Nt = N[0] * N[1] * N[2];

      data = (short***)calloc(N[0], sizeof(short**));
      for(int nx=0; nx < N[0]; nx++) {
        data[nx] = (short**)calloc(N[1], sizeof(short*));              
        if(!data[nx]) {
          cout << "failed allocation\n";
          exit(EXIT_FAILURE);
        }
        for(int ny=0; ny < N[1]; ny++) {
          data[nx][ny] = (short*)calloc(N[2], sizeof(short));              
          if(!data[nx][ny]) {
            cout << "failed allocation\n";
            exit(EXIT_FAILURE);
          }
        }
      }


      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          for(int nz=0; nz < N[2]; nz++) {
            fd.read((char*)&data[nx][ny][nz], sizeof(short));
          }
        }
      }

    }

    virtual ~ExclusionMap() {
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
      if(G[0] < 0 or G[0] >= N[0]) {
        return false;
      }
      G[1] = (int)floor(((R->y-origin[1])/delta) + 0.5);
      if(G[1] < 0 or G[1] >= N[1]) {
        return false;
      }
      G[2] = (int)floor(((R->z-origin[2])/delta) + 0.5);
      if(G[2] < 0 or G[2] >= N[2]) {
        return false;
      }

      return true;
    }

    bool onGrid(vertex *R) {
      int grid[3];
      return coordinateToGrid(R, (int*)&grid);
    }

    short value(vertex *R) {
      int grid[3];
      bool onGrid = coordinateToGrid(R, (int*)&grid);
      if(!onGrid) return 0;
      return data[grid[0]][grid[1]][grid[2]];
    }

    void translate(double dx, double dy, double dz) {
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
    }

};

#endif
