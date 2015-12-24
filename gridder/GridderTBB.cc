#include "../Main.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include "Strings.h"
#include "GAFF.h"
#include "Timer.h"

struct mytask {
  mytask(int ix, int iy, int iz, int ny, int nz, ReceptorPDBQT *rec, GAFFParameters *parm, vector< double** > *data, vector<int> *types, double padding_sqr) :_x(ix), _y(iy), _z(iz), _ny(ny), _nz(nz), _rec(rec), _parm(parm), _data(data), _types(types), _p2(padding_sqr) {
  }

  void operator()() {
    for(int i=0; i < _rec->coordinates.size(); i++) {
      vertex Rrec = _rec->coordinates[i];
      vertex dr;
      dr.x = _x - Rrec.x; 
      dr.y = _y - Rrec.y; 
      dr.z = _z - Rrec.z; 
      double dist_sqr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
      if(dist_sqr > _p2) continue;
      double dist = sqrt(dist_sqr);
      // Iterate of each ligand atom type and calc VDW/HB
      int rec_type = _rec->types[i];
      double dist_6 = dist_sqr * dist_sqr * dist_sqr;
      double dist_12 = dist_6 * dist_6;
      for(int at=0; at < _types->size(); at++) {
        int lig_type = _types->at(at);
        pair_parameter parm;
        if(lig_type > rec_type) parm = _parm->lj_map[lig_type][rec_type];
        else parm = _parm->lj_map[rec_type][lig_type];
        double du_t = (parm.A / dist_12) - (parm.B / dist_6);
        _data->at(at)[_ny][_nz] += du_t;
      }
    }
  }

  int _x, _y, _z, _ny, _nz;
  double _p2;
  ReceptorPDBQT *_rec;
  GAFFParameters *_parm;
  vector< double** > *_data;
  vector<int> *_types;
};

struct executor
{
  executor(std::vector<mytask>& t)
    :_tasks(t)
  {}
  executor(executor& e,tbb::split)
    :_tasks(e._tasks)
  {}

  void operator()(const tbb::blocked_range<size_t>& r) const {
    for (size_t i=r.begin();i!=r.end();++i)
      _tasks[i]();
  }

  std::vector<mytask>& _tasks;
};

bool getInputWithFlag(int argc, char **argv, char flag, string *value) {
  int  opt;                                                                                                                                                                                                     
  bool bopt = false;
  char gopt[2] = { flag, ':' };

  for(int i=1; i < argc; i++) {
    if(argv[i][0] == '-' && argv[i][1] == flag) {
      if(argv[i+1][0] != '-') {
        *value = argv[i+1];
        bopt = true;
        i++;
        break;
      }
    }
  }

  return bopt;
}


bool coordinateToGrid(double x, double y, double z, int *Gx, int *Gy, int *Gz, vertex *origin, double delta) {
  *Gx = (int)floor(((x-origin->x)/delta) + 0.5);
  *Gy = (int)floor(((y-origin->y)/delta) + 0.5);
  *Gz = (int)floor(((z-origin->z)/delta) + 0.5);
  return true;
}


void usage() {
  printf("Usage: gridder -d [AD4.1.DAT] -n NTHREADS(=max) -r [Receptor.PDBQT] -l [Ligand.PDBQT] (Optional: -p GRID_PADDING(=40A) -s GRID_SPACING(=0.375A) -w OutputFilenamePrefix)\n");
}


int main(int argc, char **argv) {
  string datfn, ligfn, recfn, fldfn, stoken, token, name_prefix;
  double T = 298.;
  double diel_rec = 78.5;
  double grid_resolution = 0.375;
  double padding = 40., padding_sqr;

  // GAFF parameter file
  if(!getInputWithFlag(argc, argv, 'd', &datfn)) { usage(); return -1; }
  // receptor pdbqt
  if(!getInputWithFlag(argc, argv, 'r', &recfn)) { usage(); return -1; }
  // ligand pdbqt
  if(!getInputWithFlag(argc, argv, 'l', &ligfn)) { usage(); return -1; }
  // number of processors/threads
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
    tbb::task_scheduler_init init(stringToInt(stoken));
  }
  // grid padding
  if(getInputWithFlag(argc, argv, 'p', &stoken)) {
    padding = stringToDouble(stoken);
  }
  padding_sqr = padding * padding;
  // grid spacing
  if(getInputWithFlag(argc, argv, 's', &stoken)) {
    grid_resolution = stringToDouble(stoken);
  }
  // grid padding
  if(getInputWithFlag(argc, argv, 'w', &name_prefix)) {
    cout << "> Output filename prefix: " << name_prefix << endl;
  }

  // Load AD parameters
  GAFFParameters *adp = new GAFFParameters(datfn);
  LigandPDBQT *lig = new LigandPDBQT(ligfn);

  cout << "> Types in ligand:";
  set<string>::iterator it;
  for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < adp->types.size(); j++) {
      if(adp->types[j] == *it) found = true;
    }
    if(found) cout << ' ' << *it;
    else {
      cout << "! Error: Ligand atom type " << *it << " not found in the GAFF parameters." << endl;
      exit(-1);
    }
  }
  cout << endl;

  // Load receptor file
  cout << "> Loading receptor PDBQT..." << endl;
  ReceptorPDBQT *rec = new ReceptorPDBQT(recfn, adp);
  cout << "> Receptor center: " << rec->center.x << ", " << rec->center.y << ", " << rec->center.z << endl;
  cout << "> Receptor minimum: " << rec->min.x << ", " << rec->min.y << ", " << rec->min.z << endl;
  cout << "> Receptor maximum coordinates: " << rec->max.x << ", " << rec->max.y << ", " << rec->max.z << endl;
  cout << "> Types in receptor:";
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    cout << ' ' << *it;
  }
  cout << endl;

  // Check to make sure all receptor atom types are in our parameter data set
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < adp->types.size(); j++) {
      if(adp->types[j] == *it) found = true;
    }
    if(!found) {
      cout << "! Error: Ligand atom type " << *it << " not found in the GAFF parameters." << endl;
      return EXIT_FAILURE;
    }
  }

  // Calculate grid geometries
  vertex origin = { rec->min.x - padding, rec->min.y - padding, rec->min.z - padding };
  vertex dimensions = { (rec->max.x + padding) - origin.x, (rec->max.y + padding) - origin.y, (rec->max.z + padding) - origin.z };
  int Npoints[3] = { dimensions.x / grid_resolution, dimensions.y / grid_resolution, dimensions.z / grid_resolution };
  int Ntotal = Npoints[0] * Npoints[1] * Npoints[2];

  //// Data for all grids
  // Single point calculation storage
  vector<double> u_t(lig->types_set.size());
  // Exclusion grid
  cout << "Allocating data" << endl;
  bool ***ex = (bool***)calloc(Npoints[0], sizeof(bool**));
  for(int i=0; i < Npoints[0]; i++) {
    ex[i] = (bool**)calloc(Npoints[1], sizeof(bool*));
    for(int j=0; j < Npoints[1]; j++) {
      ex[i][j] = (bool*)calloc(Npoints[2], sizeof(bool));
    }
  }

  // Dynamic storage for VDW/HB calcs
  vector<double**> data_t;
  vector<int> type_t;
  vector<ofstream*> bpm_t;
  string typefn;
  for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
    double **new_data_t = (double**)calloc(Npoints[1], sizeof(double*));
    for(int yy=0; yy < Npoints[1]; yy++) new_data_t[yy] = (double*)calloc(Npoints[2], sizeof(double));
    data_t.push_back(new_data_t);
    type_t.push_back(adp->index_for_type(*it));
    typefn = name_prefix;
    typefn.append(*it);
    typefn.append(".bpm");
    bpm_t.push_back(new ofstream(typefn, ios::out | ios::binary));
    // write header
    bpm_t[bpm_t.size()-1]->write((char*)&origin, sizeof(vertex));
    bpm_t[bpm_t.size()-1]->write((char*)&Npoints, sizeof(int)*3);
    bpm_t[bpm_t.size()-1]->write((char*)&grid_resolution, sizeof(double));
  }

  // Exclusion grid
  Timer *timer = new Timer();
  timer->start();
  cout << "> Creating exclusion grid...";
  cout.flush();
  cilk_for(int i=0; i < rec->coordinates.size(); i++) {
    vertex Rrec = rec->coordinates[i];
    double radius = rec->radii[i];
    int Grec[3];
    if(coordinateToGrid(Rrec.x, Rrec.y, Rrec.z, &Grec[0], &Grec[1], &Grec[2], &origin, grid_resolution)) {
      int searchr = ((padding + radius) / grid_resolution) + 1;
      for(int gx=(Grec[0]-searchr); gx <= Grec[0]+searchr; gx++) {
        if(gx < 0 or gx >= Npoints[0]) continue;
        for(int gy=(Grec[1]-searchr); gy <= Grec[1]+searchr; gy++) {
          if(gy < 0 or gy >= Npoints[1]) continue;
          for(int gz=(Grec[2]-searchr); gz <= Grec[2]+searchr; gz++) {
            if(gz < 0 or gz >= Npoints[2]) continue;

            double dx = (origin.x + (gx * grid_resolution)) - Rrec.x;
            double dy = (origin.y + (gy * grid_resolution)) - Rrec.y;
            double dz = (origin.z + (gz * grid_resolution)) - Rrec.z;
            if(sqrt(dx*dx + dy*dy + dz*dz) <= radius + padding) ex[gx][gy][gz] = true;
          }
        }
      }
    }
  }
  timer->stop();
  cout << "> Elapsed time (s): ";
  timer->print(&cout);
  cout << "> Starting grid calculation." << endl;

  std::vector<mytask> tasks;
  // Start a timer and start our calculations
  timer->start();
  for(int nx=0, nt=0; nx < Npoints[0]; nx++) {
    double X = (nx * grid_resolution) + origin.x;

    tasks.clear();

    for(int ny=0; ny < Npoints[1]; ny++) {
      double Y = (ny * grid_resolution) + origin.y;
      for(int nz=0; nz < Npoints[2]; nz++) {
        double Z = (nz * grid_resolution) + origin.z;

        for(int at=0; at < type_t.size(); at++)
          data_t[at][ny][nz] = 0;

        if(ex[nx][ny][nz] == false) continue;

        tasks.push_back(mytask(X, Y, Z, ny, nz, rec, adp, &data_t, &type_t, padding_sqr));
      }
    }

    executor exec(tasks);
    tbb::parallel_for(tbb::blocked_range<size_t>(0,tasks.size()), exec);

    for(int ny=0; ny < Npoints[1]; ny++) {
      for(int at=0; at < type_t.size(); at++) {
        bpm_t[at]->write((char*)data_t[at][ny], sizeof(double) * Npoints[2]);
      }
    }

    double percent_complete = nx / Npoints[0];
    cout << "\r> " << (100.*percent_complete) << " percent complete.\t\t\t" << flush;
  }
  timer->stop();
  cout << "> Elapsed time (s): ";
  timer->print(&cout);

  for(int at=0; at < type_t.size(); at++)
    bpm_t[at]->close();

  cout << endl << "> Done." << endl;
  return EXIT_SUCCESS;
}
