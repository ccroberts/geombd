
#include "Main.h"
#include "Strings.h"
#include "MartiniNB.h"
#define triMtx(i, j) (i < j ? i+j*(j-1)/2 : j+i*(i-1)/2)

enum CALCULATION_TYPE {
  CALCULATION_ES,
  CALCULATION_VDW
};


int main(int argc, char **argv) {
  if(argc < 7) {
    errorf("Usage: Gridder ES  cgPQR dxOUT GRIDRESOLUTION THREADS PADDING TEMP SOLVENTDIELECTRIC CONC+1 CONC-1 CONC+2 CONC-2\n       Gridder VDW cgPQR dxOUT GRIDRESOLUTION THREADS PADDING PROBE_TYPE\n");
  }

  __cilkrts_set_param("nworkers", argv[5]);
  int numWorkers = __cilkrts_get_nworkers();
  printf("> Number of threads: %d\n", numWorkers);

  // collect arguments
  CALCULATION_TYPE calculationType = (argv[1][0] == 'E') ? CALCULATION_ES : CALCULATION_VDW;
  string filenameCG(argv[2]);
  string filenameDX(argv[3]);
  double resolution = charToDouble(argv[4]);
  double pad = charToDouble(argv[6]);
  double pad2 = pad*pad;
  string filenameENV("");

  int vdwProbeType = 12;
  double T = 298.;
  double epsR = 78.5;
  double ConcMonoCat = 0.1;//M
  double ConcMonoAni = 0.1;//M
  double ConcDiCat = 0.0;//M
  double ConcDiAni = 0.0;//M
  if(calculationType == CALCULATION_VDW) {
    vdwProbeType = atoi(argv[7]);
    if(argc > 9) {
      filenameENV = argv[8];
    }
  } else {
    T = charToDouble(argv[7]);
    epsR = charToDouble(argv[8]);
    ConcMonoCat = charToDouble(argv[9]);
    ConcMonoAni = charToDouble(argv[10]);
    ConcDiCat   = charToDouble(argv[11]);
    ConcDiAni   = charToDouble(argv[12]);
    if(argc > 13) {
      filenameENV = argv[13];
    }
  }

  // calculate atomic system parameters
  double X1 =  1e10;
  double X2 = -1e10;
  double Y1 =  1e10;
  double Y2 = -1e10;
  double Z1 =  1e10;
  double Z2 = -1e10;

  vector<double> aX, aY, aZ, aQ;
  int Natoms = 0;
  string line, token;
  double dtoken;
  ifstream fd(filenameCG.c_str(), ifstream::in);
  while(getline(fd, line)) {
    parseNextValue(&line, &token);
    if(token == "ATOM") {
      line = line.substr(24, line.length() - 24);
      parseNextValue(&line, &token);
      dtoken = stringToDouble(token);
      aX.push_back(dtoken);
      if(dtoken < X1) X1 = dtoken;
      if(dtoken > X2) X2 = dtoken;
      parseNextValue(&line, &token);
      dtoken = stringToDouble(token);
      aY.push_back(dtoken);
      if(dtoken < Y1) Y1 = dtoken;
      if(dtoken > Y2) Y2 = dtoken;
      parseNextValue(&line, &token);
      dtoken = stringToDouble(token);
      aZ.push_back(dtoken);
      if(dtoken < Z1) Z1 = dtoken;
      if(dtoken > Z2) Z2 = dtoken;
      parseNextValue(&line, &token);
      if(calculationType == CALCULATION_ES)
        aQ.push_back(stringToDouble(token)); //charge
      parseNextValue(&line, &token);
      if(calculationType == CALCULATION_VDW)
        aQ.push_back(stringToDouble(token)); //lj bead type
      Natoms ++;
    }
  }
  fd.close();

  if(filenameENV != "") {
    cout << "> Loading in environmental atoms" << endl;
    double EX, EY, EZ;
    ifstream fd(filenameENV.c_str(), ifstream::in);
    while(getline(fd, line)) {
      parseNextValue(&line, &token);
      if(token == "ATOM") {
        line = line.substr(24, line.length() - 24);
        parseNextValue(&line, &token);
        EX = stringToDouble(token);
        if(EX < X1 - pad) continue;
        if(EX > X2 + pad) continue;
        parseNextValue(&line, &token);
        EY = stringToDouble(token);
        if(EY < Y1 - pad) continue;
        if(EY > Y2 + pad) continue;
        parseNextValue(&line, &token);
        EZ = stringToDouble(token);
        if(EZ < Z1 - pad) continue;
        if(EZ > Z2 + pad) continue;
        aX.push_back(EX);
        aY.push_back(EY);
        aZ.push_back(EZ);
        parseNextValue(&line, &token);
        if(calculationType == CALCULATION_ES)
          aQ.push_back(stringToDouble(token)); //charge
        parseNextValue(&line, &token);
        if(calculationType == CALCULATION_VDW)
          aQ.push_back(stringToDouble(token)); //lj bead type
        Natoms ++;
      }
    }
    fd.close();
  }

  printf("> %d atoms in molecular system.\n", Natoms);

  // pad bounds
  X1 -= pad;
  Y1 -= pad;
  Z1 -= pad;
  X2 += pad;
  Y2 += pad;
  Z2 += pad;

  // screening constant
  double ions = (/*1+*/1. * (ConcMonoCat * 6.022e23 / 1e27))
              + (/*1-*/1. * (ConcMonoAni * 6.022e23 / 1e27))
              + (/*2+*/4. * (ConcDiCat * 6.022e23 / 1e27))
              + (/*2-*/4. * (ConcDiAni * 6.022e23 / 1e27));
  double kappa = sqrt((4. * M_PI * kC * ions) / (epsR * kB * T));
  if(calculationType == CALCULATION_ES) {
    printf("> Screening parameter, Kappa: %f\n", kappa);
  }


  // calculate grid parameters
  int NpointsX = (X2-X1) / resolution;
  int NpointsY = (Y2-Y1) / resolution;
  int NpointsZ = (Z2-Z1) / resolution;
  int Npoints = NpointsX * NpointsY * NpointsZ;
  printf("> Grid points( X: %d Y: %d Z: %d Total: %d )\n", NpointsX, NpointsY, NpointsZ, Npoints);
  printf("> Estimated memory size: %fGB\n", Npoints*sizeof(double)*1e-9);


  // allocate memory
  double ***data = (double***)calloc(NpointsX, sizeof(double**));
  if(data == NULL) errorf("Error allocating memory.\n");
  for(int nx=0; nx < NpointsX; nx++) {
    data[nx] = (double**)calloc(NpointsY, sizeof(double*));              
    if(data[nx] == NULL) errorf("Error allocating memory.\n");
    for(int ny=0; ny < NpointsY; ny++) {
      data[nx][ny] = (double*)calloc(NpointsZ, sizeof(double));              
      if(data[nx][ny] == NULL) errorf("Error allocating memory.\n");
    }
  }


  // calculate potential across grid
  if(calculationType == CALCULATION_ES)
    printf("> Starting electrostatic grid potential calculation.\n");
  else
    printf("> Starting VDW grid potential calculation.\n");


  for(int nx=0; nx < NpointsX; nx++) {
    double X = (nx * resolution) + X1;
    for(int ny=0; ny < NpointsY; ny++) {
      double Y = (ny * resolution) + Y1;
      cilk_for(int nz=0; nz < NpointsZ; nz++) {
        double Z = (nz * resolution) + Z1;

        double dr[3], u = 0.;
        for(int i=0; i < aZ.size(); i++) {
          double du = 0.;

          if(calculationType == CALCULATION_ES) {
            if(aQ[i] == 0.) continue;
            dr[0] = X - aX[i];
            if(dr[0] > pad) continue;
            dr[1] = Y - aY[i];
            if(dr[1] > pad) continue;
            dr[2] = Z - aZ[i];
            if(dr[2] > pad) continue;
            double distSqr = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if(distSqr > pad2) continue;
            if(distSqr < 2.25) {
              du = (aQ[i] > 0.) ? 1.5e3 : -1.5e3;
              break;
            } else {
              double dist = sqrt(distSqr);
              du = kC * aQ[i] * exp(-dist * kappa) / (epsR * dist);
              if(du > 1.5e3) du = 1.5e3;
              if(du < -1.5e3) du = -1.5e3;
            }
          } else {
            dr[0] = X - aX[i];
            if(dr[0] > pad) continue;
            dr[1] = Y - aY[i];
            if(dr[1] > pad) continue;
            dr[2] = Z - aZ[i];
            if(dr[2] > pad) continue;
            double distSqr = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if(distSqr > pad2) continue;
            if(distSqr < 2.25) {
              du = 1.5e3;
              break;
            } else {
              double invDist2 = 1./distSqr;
              double invDist6 = invDist2 * invDist2 * invDist2;

              du = (c12[triMtx(vdwProbeType, (int)aQ[i])] * invDist6 * invDist6 - c6[triMtx(vdwProbeType, (int)aQ[i])] * invDist6);
              if(distSqr >= 81.) {
                double dist = sqrt(distSqr);
                du *= (1 - pow((dist - 9.) / 3., 2));
              }
              if(du > 1.5e3) du = 1.5e3;
              if(du < -1.5e3) du = -1.5e3;
            }
          }
          u += du;
        }
        data[nx][ny][nz] += u;
      }
    }
    printf("> %f%% complete.\n", 100. * (((double)nx+1) / (double)NpointsX));
  }

  printf("> Starting to write OpenDX potential map...\n");
  // output opendx data
  FILE *fdo;
  fdo = fopen(filenameDX.c_str(), "w");
  fprintf(fdo, "object 1 class gridpositions counts %d %d %d\n", NpointsX, NpointsY, NpointsZ);
  fprintf(fdo, "origin %12.6e %12.6e %12.6e\n", X1, Y1, Z1);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", resolution, 0., 0.);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., resolution, 0.);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., 0., resolution);
  fprintf(fdo, "object 2 class gridconnections counts %d %d %d\n", NpointsX, NpointsY, NpointsZ);
  fprintf(fdo, "object 3 class array type double rank 0 items %d data follows\n", Npoints);
  int i[3] = { 0, 0, 0 };
  for(int it=0; it < Npoints; it++) {
    fprintf(fdo, "%12.6e ", data[i[0]][i[1]][i[2]]);
    if((it+1) % 3 == 0) fprintf(fdo, "\n");
    i[2]++;
    if(i[2] >= NpointsZ) { i[2] = 0; i[1]++; }
    if(i[1] >= NpointsY) { i[1] = 0; i[0]++; }
  }
  if(Npoints % 3 != 0) fprintf(fdo, "\n");
  fprintf(fdo, "attribute \"dep\" string \"positions\"\n");
  fprintf(fdo, "object \"regular positions regular connections\" class field\n");
  fprintf(fdo, "component \"positions\" value 1\n");
  fprintf(fdo, "component \"connections\" value 2\n");
  fprintf(fdo, "component \"data\" value 3\n");
  fclose(fdo);

  // free memory
  for(int nx=0; nx < NpointsX; nx++) {
    for(int ny=0; ny < NpointsY; ny++) {
      free(data[nx][ny]);
    }
    free(data[nx]);
  }
  free(data);

  printf("> Done.\n");
  return EXIT_SUCCESS;
}
