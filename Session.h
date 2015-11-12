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
    int Nbind;
    int Nexit;
    int Ntlim;
    double Davg;
    vector< Body* > conformations;
    vector< Body* > ligands;
    vector< BindingCriteria* > bindingCriteria;

  public:
    Session(Model *m, SimulationConfig s);
    ~Session();

    virtual void populateLigands();
    virtual void printRateConstant();

};


class SessionRadial : public Session {
  friend class Model;
  protected:
    double b;
    double q, q2;

  public:
    SessionRadial(Model *m);

    virtual void populateLigands();
    virtual void printRateConstant();

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

    virtual void populateLigands();
    virtual void printRateConstant();

};


class SessionAbsoluteRadial : public Session {
  friend class Model;
  protected:
    vertex start;
    double b;
    double q, q2;

  public:
    SessionAbsoluteRadial(Model *m);

    virtual void populateLigands();
    virtual void printRateConstant();

};


#endif
