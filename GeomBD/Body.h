#ifndef Body_h_
#define Body_h_

#include "Main.h"
#include "Bead.h"

 
class Model;
class Session;


 
class Body {
  public:
    Model *model;
    Session *session;

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

  public:
    double t_dwell, t_dwell_max, t_dwell_total;

  public:
    Body();
    Body(Model *m, Session *s);
    virtual ~Body();

  public:
    virtual void define();
    void writePDB(fstream &outf, char chain);

  public:
    virtual void translate(double dx, double dy, double dz, bool suppressWarning=false);
    virtual void center();
    virtual void rotate(double dax, double day, double daz);

};


#endif
