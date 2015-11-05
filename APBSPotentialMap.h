#ifndef _APBSPotentialMap_h_
#define _APBSPotentialMap_h_

#include "Main.h"
#include "Strings.h"

 
class APBSPotentialMap {
  protected:
    bool header = true,
         body = false,
         footer = false;

    int N[3];
    int Nt;
    int i[3];
    int it;
    double origin[3];
    double delta[3];
    double ***data;

  public:
    APBSPotentialMap(string dx_filename, double KbT) {
      ifstream fd;
      string line, token;

      for(int c=0; c < 3; c++) i[c] = 0.;
      it = 0.;

      fd.open(dx_filename.c_str(), ifstream::in);
      while(getline(fd, line)) {
        if(header) {
          parseNextValue(&line, &token);
          if(token == "object") {
            parseNextValue(&line, &token);
            if(token == "1") {
              for(int c=0; c < 3; c++)
                parseNextValue(&line, &token);
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
            if(token == "3") {
              header = false;
              body = true;
            }
          }
          if(token == "origin") {
            for(int c=0; c < 3; c++) {
              parseNextValue(&line, &token);
              origin[c] = stringToDouble(token);
            }
          }
          if(token == "delta") {
            double d[3];
            for(int c=0; c < 3; c++) {
              parseNextValue(&line, &token);
              d[c] = stringToDouble(token);
              if(d[c] != 0.) delta[c] = d[c];
            }
          }
          continue;
        }
        if(body) {
          while(parseNextValue(&line, &token)) {
            double val = stringToDouble(token) * KbT;
            data[i[0]][i[1]][i[2]] = val;
            i[2]++;
            it++;
            if(i[2] >= N[2]) { i[2] = 0; i[1]++; }
            if(i[1] >= N[1]) { i[1] = 0; i[0]++; }
            if(it >= Nt) { body = false; footer = true; }
          }

          continue;
        }
      }
    }

    virtual ~APBSPotentialMap() {
      for(int nx=0; nx < N[0]; nx++) {
        for(int ny=0; ny < N[1]; ny++) {
          free(data[nx][ny]);
        }
        free(data[nx]);
      }
      free(data);
    }

    bool coordinateToGrid(vertex *R, int *G) {
      G[0] = (int)floor(((R->x-origin[0])/delta[0]) + 0.5);
      if(G[0] <= 0 or G[0] >= N[0]-1) {
        return false;
      }
      G[1] = (int)floor(((R->y-origin[1])/delta[1]) + 0.5);
      if(G[1] <= 0 or G[1] >= N[1]-1) {
        return false;
      }
      G[2] = (int)floor(((R->z-origin[2])/delta[2]) + 0.5);
      if(G[2] <= 0 or G[2] >= N[2]-1) {
        return false;
      }

      return true;
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
      dU[0] = (data[grid[0]-1][grid[1]][grid[2]] - data[grid[0]+1][grid[1]][grid[2]]) / (2. * delta[0]);
      dU[1] = (data[grid[0]][grid[1]-1][grid[2]] - data[grid[0]][grid[1]+1][grid[2]]) / (2. * delta[1]);
      dU[2] = (data[grid[0]][grid[1]][grid[2]-1] - data[grid[0]][grid[1]][grid[2]+1]) / (2. * delta[2]);

      if(e) *e = E * q;

      F->x += q * dU[0];
      F->y += q * dU[1];
      F->z += q * dU[2];

      return true;
    }

    void translate(double dx, double dy, double dz) {
      origin[0] += dx;
      origin[1] += dy;
      origin[2] += dz;
    }

};

#endif
