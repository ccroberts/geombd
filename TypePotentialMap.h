#ifndef _TypePotentialMap_h_
#define _TypePotentialMap_h_

#include "Main.h"
#include "Strings.h"
#include "PotentialMap.h"

 
class TypePotentialMap : public PotentialMap {
  public:
    char type;

  public:
    TypePotentialMap(string dx_filename, char atomtype) : PotentialMap(dx_filename) {
      type = atomtype;
    }

};

#endif
