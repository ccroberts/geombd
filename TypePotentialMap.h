#ifndef _TypePotentialMap_h_
#define _TypePotentialMap_h_

#include "Main.h"
#include "Strings.h"
#include "BinaryPotentialMap.h"

 
class TypePotentialMap : public BinaryPotentialMap {
  public:
    TypePotentialMap(string bpm_filename, string atomtype) : BinaryPotentialMap(bpm_filename, atomtype) {
    }

};

#endif
