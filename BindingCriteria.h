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


  public:
    BindingCriteria() {
      Nbind.set_value(0);
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
      for(int i=0; i < pairs.size(); i++) {
        Bead *lb = ligand->beads[pairs[i].laid];
        vertex dr = { lb->R.x - pairs[i].target.x,
                      lb->R.y - pairs[i].target.y,
                      lb->R.z - pairs[i].target.z };
        double l2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        if(l2 > pairs[i].r2) return false;
      }
      return true;
    }


};


#endif
