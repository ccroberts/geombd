
#include "Model.h"
#include "Strings.h"



void writeBead(Bead *bi, unt index, char chain, fstream &outf) {
  outf << "ATOM  ";
  outf.width(5);
  outf << right << index;
  outf << "      ";
  outf.width(3);
  outf << " P " << " ";
  outf << chain;
  outf << "        ";
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.x/1.;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.y/1.;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.z/1.;
  outf << " ";
  outf.setf(ios::right, ios::adjustfield);
  outf.precision(4);
  outf.width(7);
  outf << bi->q << " " << bi->type << endl;
}



void Model::writeCoordinatesPQR() {
  // create output file
  fstream outf(ofn.c_str(), ios::out | ios::app);
  int index = 0;
  char chain = 'A';

  for(int b=0; b < ligands.size(); b++) {
    Body *ligand = ligands[b];

    ligand->writePDB(outf, chain);
    /*
    for(int i=0; i < ligand->beads.size(); i++) {
      index++;
      if(index >= 100000) index = 1;

      Bead *bi = ligand->beads[i];
      writeBead(bi, index, chain, outf);
    }
    */

    chain++;
    if(chain > 'z') chain = 'A';
  }
  outf << "END" << endl;
  outf.close();
}



