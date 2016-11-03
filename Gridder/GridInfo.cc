#include "Gridder.h"


int main(int argc, char **argv) {
  double origin[3];
  int N[3];
  double delta;

  ifstream fd(argv[1], ios::in | ios::binary);

  fd.read((char*)origin, sizeof(double) * 3);
  fd.read((char*)N, sizeof(int) * 3);
  fd.read((char*)&delta, sizeof(double));

  printf("Gridder Binary Grid Information:\n");
  printf("- Origin:       %10.4f %10.4f %10.4f\n", origin[0], origin[1], origin[2]);
  printf("- Points:       %10d %10d %10d\n", N[0], N[1], N[2]);
  printf("- Spacing:      %10.4f\n", delta);
  printf("- Dimensions:   %10.4f %10.4f %10.4f\n", N[0]*delta, N[1]*delta, N[2]*delta);
}
