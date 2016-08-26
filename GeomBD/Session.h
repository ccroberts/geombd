#ifndef Session_h_
#define Session_h_

#include "Main.h"

class Model;
class Body;
class BindingCriteria;

 
class Session {
  friend class Model;

  protected:
    int id;
    Model *model;
    SimulationConfig type;
    int Nreplicates;
    cilk::reducer< cilk::op_add<int> > Nbind;
    cilk::reducer< cilk::op_add<int> > Nexit;
    cilk::reducer< cilk::op_add<int> > Ntlim;
    cilk::reducer< cilk::op_add<double> > t_avgt;
    double Davg;
    vector< Body* > conformations;
    vector< Body* > ligands;
    vector< BindingCriteria* > bindingCriteria;
    deque<double> beta_history;
    bool done;

  public:
    Session(Model *m, SimulationConfig s);
    ~Session();

    virtual void populateLigands();
    virtual void positionLigand(Body *body) { };
    virtual void printRateConstant() { };
    virtual void checkLigand(Body *body) { };
    virtual void recordBeta(double beta);
    virtual void checkConvergence();
    virtual void finalize() { };

};


class SessionRadial : public Session {
  friend class Model;
  protected:
    double b;
    double q, q2;

  public:
    SessionRadial(Model *m);

    virtual void positionLigand(Body *body);
    virtual void printRateConstant();
    virtual void checkLigand(Body *body);

};


class SessionAbsolutePeriodic : public Session {
  friend class Model;
  protected:
    vertex start;
    vertex bounds;
    double b;
    double t_max;

  public:
    SessionAbsolutePeriodic(Model *m);

    virtual void positionLigand(Body *body);
    virtual void printRateConstant();
    virtual void checkLigand(Body *body);

};


class SessionAbsoluteRadial : public Session {
  friend class Model;
  protected:
    vertex start;
    double b;
    double q, q2;

  public:
    SessionAbsoluteRadial(Model *m);

    virtual void positionLigand(Body *body);
    virtual void printRateConstant();
    virtual void checkLigand(Body *body);

};


#endif
