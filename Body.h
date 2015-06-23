#ifndef Body_h_
#define Body_h_

#include "Main.h"
#include "Bead.h"

 
class Model;


 
class Body {
  public:
    Model *model;

  public:
    int N;
    vector< Bead* > beads;

  public:
    double m, I;
    double D, Da;
    vertex R, Ra;
    vertex F, Fa;
    double r, r_max;

  public:
    double t;
    double dt;
    bool done;     //no longer simulating
    bool bound;    //bound to binding site
    bool bulk;     //left the simulation volume

  public:
    Body();
    Body(Model *mod);
    virtual ~Body();

  public:
    virtual void define();
    void writePDB(string filename);

  public:
    virtual void translate(double dx, double dy, double dz);
    virtual void center();
    virtual void rotate(double dax, double day, double daz);

};


#endif
