#ifndef BindingCriteria_h_
#define BindingCriteria_h_
#include <cilk/reducer_opadd.h>

struct BindingPair {
  double r2;
  vertex target;
  int laid;
};


class BindingCriteria {
  public:
    vector<BindingPair> pairs;
    cilk::reducer< cilk::op_add<int> > Nbind;
    cilk::reducer< cilk::op_add<double> > t_avgt;
    bool AND;


  public:
    BindingCriteria(bool cAND=true) {
      Nbind.set_value(0);
      t_avgt.set_value(0.);
      AND = cAND;
    }

    ~BindingCriteria() {}

    void addPair(double rx, double ry, double rz, int lig_id, double r) {
      BindingPair bp;
      bp.r2 = r*r;
      bp.target.x = rx;
      bp.target.y = ry;
      bp.target.z = rz;
      bp.laid = lig_id;
      pairs.push_back(bp);
    }

    bool checkBinding(Body *ligand) {
      if(AND) return checkBindingAND(ligand);
      return checkBindingOR(ligand);
    }

    bool checkBindingAND(Body *ligand) {
      for(int i=0; i < pairs.size(); i++) {
        vertex dr;
        if(pairs[i].laid < 0) {
          dr.x = ligand->R.x - pairs[i].target.x;
          dr.y = ligand->R.y - pairs[i].target.y;
          dr.z = ligand->R.z - pairs[i].target.z;
        } else {
          Bead *lb = ligand->beads[pairs[i].laid];
          dr.x = lb->R.x - pairs[i].target.x;
          dr.y = lb->R.y - pairs[i].target.y;
          dr.z = lb->R.z - pairs[i].target.z;
        }
        double l2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        if(l2 > pairs[i].r2) return false;
      }
      return true;
    }

    bool checkBindingOR(Body *ligand) {
      for(int i=0; i < pairs.size(); i++) {
        vertex dr;
        if(pairs[i].laid < 0) {
          dr.x = ligand->R.x - pairs[i].target.x;
          dr.y = ligand->R.y - pairs[i].target.y;
          dr.z = ligand->R.z - pairs[i].target.z;
        } else {
          Bead *lb = ligand->beads[pairs[i].laid];
          dr.x = lb->R.x - pairs[i].target.x;
          dr.y = lb->R.y - pairs[i].target.y;
          dr.z = lb->R.z - pairs[i].target.z;
        }
        double l2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        if(l2 <= pairs[i].r2) return true;
      }
      return false;
    }


};


#endif
