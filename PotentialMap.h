#ifndef _PotentialMap_h_
#define _PotentialMap_h_

#include "Main.h"
#include "Strings.h"

 
class PotentialMap {
  protected:
    bool header = true,
         body = false,
         footer = false;

    int N[3];
    int Nt;
    int i[3];
    int it;
    double origin[3];
    double delta;
    double ***data;

  public:
    PotentialMap(string dx_filename) {
      ifstream fd;
      string line, token;

      for(int c=0; c < 3; c++) i[c] = 0.;
      it = 0.;

      fd.open(dx_filename.c_str(), ifstream::in);
      while(getline(fd, line)) {
        if(header) {
          parseNextValue(&line, &token);
          if(token == "NELEMENTS") {
            for(int c=0; c < 3; c++) {
              parseNextValue(&line, &token);
              N[c] = stringToInt(token);
            }
            Nt = N[0]*N[1]*N[2];
            data = (double***)calloc(N[0], sizeof(double**));
            for(int nx=0; nx < N[0]; nx++) {
              data[nx] = (double**)calloc(N[1], sizeof(double*));              
              for(int ny=0; ny < N[1]; ny++) {
                data[nx][ny] = (double*)calloc(N[2], sizeof(double));              
              }
            }
          }
          if(token == "SPACING") {
            parseNextValue(&line, &token);
            delta = stringToDouble(token);
          }
          if(token == "CENTER") {
            for(int c=0; c < 3; c++) {
              parseNextValue(&line, &token);
              origin[c] = stringToDouble(token);
              origin[c] -= (0.5 * N[c]) * delta;
            }
            header = false;
            body = true;
          }
          continue;
        }
        if(body) {
          parseNextValue(&line, &token);
          double val = stringToDouble(token);
          data[i[0]][i[1]][i[2]] = val;
          i[0]++;
          it++;
          if(i[0] >= N[0]) { i[0] = 0; i[1]++; }
          if(i[1] >= N[1]) { i[1] = 0; i[2]++; }
          if(it >= Nt) { body = false; }
          continue;
        }

      }
    }

    virtual ~PotentialMap() {
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

      double dU1[3], dU2[3];
      dU1[0] = (data[grid[0]-1][grid[1]][grid[2]] - E) / (delta);
      dU1[1] = (data[grid[0]][grid[1]-1][grid[2]] - E) / (delta);
      dU1[2] = (data[grid[0]][grid[1]][grid[2]-1] - E) / (delta);

      dU2[0] = (E - data[grid[0]+1][grid[1]][grid[2]]) / (delta);
      dU2[1] = (E - data[grid[0]][grid[1]+1][grid[2]]) / (delta);
      dU2[2] = (E - data[grid[0]][grid[1]][grid[2]+1]) / (delta);

      *e = E;

      F->x += ((dU1[0] + dU2[0]) / 2.);
      F->y += ((dU1[1] + dU2[1]) / 2.);
      F->z += ((dU1[2] + dU2[2]) / 2.);

      return true;
    }

    void translate(double dx, double dy, double dz) {
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
    }

};

#endif
