#ifndef Body_h_
#define Body_h_

#include "Main.h"
#include "Bead.h"

 
class Model;
class Session;


 
class Body {
  private:
    vertex _R, _Ra;
    vertex _F, _Fa;
    double _t;

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
    double mF;
    double r, r_max;

  public:
    double t;
    double dt;

  public:
    double t_dwell, t_dwell_max, t_dwell_total;

  public:
    Body();
    Body(Model *m, Session *s);
    virtual ~Body();

  public:
    virtual void define();
    virtual void save();
    virtual void restore();

  public:
    virtual bool translate(double dx, double dy, double dz);
    virtual void center();
    virtual void rotate(double dax, double day, double daz);

  public:
    void writePDB(fstream &outf, char chain);

};


#endif
