#ifndef Bead_h_
#define Bead_h_

#include "Main.h"


 
class Bead {
  public:
    vertex R;
    vertex F;
    double q;
    double r;
    double m;
    string type;

  public:
    Bead() {
      R.x = 0.;
      R.y = 0.;
      R.z = 0.;
      F.x = 0.;
      F.y = 0.;
      F.z = 0.;
      q = 0.;
      r = 0.;
      m = 0.;
    }


    ~Bead() {
    }


    void translate(double dx, double dy, double dz) {
      R.x += dx;
      R.y += dy;
      R.z += dz;
    }


    void rotate(double dax, double day, double daz) {
      double rm[3][3];
      double cost, sint;

      double v[3] = { R.x, R.y, R.z };
      double p[3] = { 0., 0., 0. };

      if(dax != 0) {
        cost = cos(dax);
        sint = sin(dax);

        rm[0][0] = 1;
        rm[0][1] = 0.;
        rm[0][2] = 0.;

        rm[1][0] = 0.;
        rm[1][1] = cost;
        rm[1][2] = -sint;

        rm[2][0] = 0.;
        rm[2][1] = sint;
        rm[2][2] = cost;

        for(int n=0; n < 3; n++)
          for(int m=0; m < 3; m++)
            p[n] += rm[n][m] * v[m];

        for(int n=0; n < 3; n++) {
          v[n] = p[n];
          p[n] = 0.;
        }
      }

      if(day != 0) {
        cost = cos(day);
        sint = sin(day);

        rm[0][0] = cost;
        rm[0][1] = 0.;
        rm[0][2] = sint;

        rm[1][0] = 0.;
        rm[1][1] = 1;
        rm[1][2] = 0.;

        rm[2][0] = -sint;
        rm[2][1] = 0.;
        rm[2][2] = cost;

        for(int n=0; n < 3; n++)
          for(int m=0; m < 3; m++)
            p[n] += rm[n][m] * v[m];

        for(int n=0; n < 3; n++) {
          v[n] = p[n];
          p[n] = 0.;
        }
      }

      if(daz != 0) {
        cost = cos(daz);
        sint = sin(daz);

        rm[0][0] = cost;
        rm[0][1] = -sint;
        rm[0][2] = 0.;

        rm[1][0] = sint;
        rm[1][1] = cost;
        rm[1][2] = 0.;

        rm[2][0] = 0.;
        rm[2][1] = 0.;
        rm[2][2] = 1;

        for(int n=0; n < 3; n++)
          for(int m=0; m < 3; m++)
            p[n] += rm[n][m] * v[m];

        for(int n=0; n < 3; n++) {
          v[n] = p[n];
          p[n] = 0.;
        }
      }

      R.x = v[0];
      R.y = v[1];
      R.z = v[2];
    }


    void rotateAbout(double ox, double oy, double oz, double dax, double day, double daz) {
      translate(-ox, -oy, -oz);
      rotate(dax, day, daz);
      translate(ox, oy, oz);
    }

};


#endif
