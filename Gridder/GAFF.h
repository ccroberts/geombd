#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;


struct pair_parameter {
  double A;
  double B;
  double e;
  double s;
};

typedef vector<pair_parameter> Vec;
typedef vector<Vec> PairMap;


class GAFFParameters {
  public:
    /*
    double w_vdw;
    double w_hbond;
    double w_elec;
    double w_desolv;
    */

    vector<string> types;
    vector<double> Rii;
    vector<double> epsii;
    /*
    vector<double> vol;
    vector<double> solpar;
    vector<double> Rij_hb;
    vector<double> epsij_hb;
    vector<int> hbond_type;
    */

    PairMap lj_map;

  public:
    GAFFParameters(string filename) {
      string line, token;
      ifstream fd;

      fd.open(filename.c_str(), ifstream::in);
      while(getline(fd, line)) {
        if(line[0] == '#') continue;
        parseNextValue(&line, &token);

        if(token.length() > 2) {
          // free energy weights
          /*
          if(token.substr(0, 2) == "FE") {
            char weight_type = token.substr(9, 1)[0];
            parseNextValue(&line, &token);
            double value = stringToDouble(token);
            switch(weight_type) {
              case 'v':
                w_vdw = value;
                break;
              case 'h':
                w_hbond = value;
                break;
              case 'e':
                w_elec = value;
                break;
              case 'd':
                w_desolv = value;
                break;
            }
          }
          */
          // everything else
          if(token == "atom_par") {
            parseNextValue(&line, &token);
            types.push_back(token);
            parseNextValue(&line, &token);
            Rii.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            epsii.push_back(stringToDouble(token));
            /*
            parseNextValue(&line, &token);
            vol.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            solpar.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            Rij_hb.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            epsij_hb.push_back(stringToDouble(token));
            parseNextValue(&line, &token);
            hbond_type.push_back(stringToInt(token));
            */
          }
        }
      }

      lj_map = PairMap(types.size());

      for(int i=0; i < types.size(); i++) {
        for(int j=0; j < i+1; j++) {
          pair_parameter parm;

          double rij = 0.5 * (Rii[i] + Rii[j]);
          double eij = sqrt(epsii[i] * epsii[j]);
          parm = {
            /*A*/4. * eij * pow(rij, 12),
            /*B*/4. * eij * pow(rij, 6),
            eij,
            rij
          };

          lj_map[i].push_back(parm);
        }
      }
    }

    int index_for_type(string type) {
      for(int i=0; i < types.size(); i++) {
        if(types[i] == type) return i;
      }
      cout << "* Warning: GAFFParameters::index_for_type returning -1" << endl;
      return -1;
    }


};


class LigandPDBQT {
  public:
    set<string> types_set;

  public:
    LigandPDBQT(string filename) {
      string line, token;
      ifstream fd;

      fd.open(filename.c_str(), ifstream::in);
      while(getline(fd, line)) {
        if(line[0] == '#') continue;
        if(line.substr(0, 4) == "ATOM") {
          string t = line.substr(78, 2);
          string tt = rtrim(t);
          if(tt == "A") tt = "C";
          if(tt == "HD") tt = "H";
          if(tt == "HS") tt = "H";
          if(tt == "NA") tt = "N";
          if(tt == "NS") tt = "N";
          if(tt == "OA") tt = "O";
          if(tt == "OS") tt = "O";
          if(tt == "SA") tt = "S";
          if(tt == "CL") tt = "Cl";
          if(tt == "BR") tt = "Br";
          if(tt == "MG") tt = "Mg";
          if(tt == "CA") tt = "Ca";
          if(tt == "MN") tt = "Mn";
          if(tt == "FE") tt = "Fe";
          if(tt == "ZN") tt = "Zn";
          types_set.insert(tt);
        }
      }
    }
};



class ReceptorPDBQT {
  public:
    vector<vertex> coordinates;
    vertex center;
    vertex max;
    vertex min;
    vector<int> types;
    set<string> types_set;
    vector<double> charges;
    vector<double> radii;

  public:
    ReceptorPDBQT(string filename, GAFFParameters *adp) {
      string line, token;
      ifstream fd;

      center.x = center.y = center.z = 0.;
      max.x = max.y = max.z = -1e9;
      min.x = min.y = min.z =  1e9;

      fd.open(filename.c_str(), ifstream::in);
      while(getline(fd, line)) {
        if(line[0] == '#') continue;
        if(line.substr(0, 4) == "ATOM" or line.substr(0, 6) == "HETATM") {
          // coordinates
          vertex R;
          R.x = stringToDouble(line.substr(30, 8));
          R.y = stringToDouble(line.substr(38, 8));
          R.z = stringToDouble(line.substr(46, 8));
          coordinates.push_back(R);
          center.x += R.x;
          center.y += R.y;
          center.z += R.z;
          if(R.x > max.x) max.x = R.x;
          if(R.y > max.y) max.y = R.y;
          if(R.z > max.z) max.z = R.z;
          if(R.x < min.x) min.x = R.x;
          if(R.y < min.y) min.y = R.y;
          if(R.z < min.z) min.z = R.z;
          // type
          string t = line.substr(78, 2);
          string tt = rtrim(t);
          if(tt == "A") tt = "C";
          if(tt == "HD") tt = "H";
          if(tt == "HS") tt = "H";
          if(tt == "NA") tt = "N";
          if(tt == "NS") tt = "N";
          if(tt == "OA") tt = "O";
          if(tt == "OS") tt = "O";
          if(tt == "SA") tt = "S";
          if(tt == "CL") tt = "Cl";
          if(tt == "BR") tt = "Br";
          if(tt == "MG") tt = "Mg";
          if(tt == "CA") tt = "Ca";
          if(tt == "MN") tt = "Mn";
          if(tt == "FE") tt = "Fe";
          if(tt == "ZN") tt = "Zn";
          if(adp)
            types.push_back(adp->index_for_type(tt));
          types_set.insert(tt);
          // charge
          double q = stringToDouble(line.substr(69, 8));
          charges.push_back(q);
          // radius
          int ati = -1;
          double radius = 1.4;
          if(adp) {
            ati = adp->index_for_type(tt);
            if(ati >= 0) {
              radius = adp->Rii[ati] / 2.;
            }
          }
          radii.push_back(radius);
        }
      }

      center.x /= coordinates.size();
      center.y /= coordinates.size();
      center.z /= coordinates.size();
    }
};
