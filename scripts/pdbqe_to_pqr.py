#!/usr/bin/env python3

import sys


def run(params, pdbqefn):
  t2r = {}

  for line in open(sys.argv[1], 'r'):
    if line.startswith('atom_par'):
      sp = line.split()
      t = sp[1]
      r = float(sp[2]) / 2.
      t2r[t] = r

  for line in open(pdbqefn, 'r'):
    if line.startswith('ATOM') or line.startswith('HETATM'):
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      q = float(line[70:76])
      tt = line[77:].strip()
      if tt == "A": tt = "C";
      if tt == "HD": tt = "H";
      if tt == "HS": tt = "H";
      if tt == "NA": tt = "N";
      if tt == "NS": tt = "N";
      if tt == "OA": tt = "O";
      if tt == "OS": tt = "O";
      if tt == "SA": tt = "S";
      if tt == "CL": tt = "Cl";
      if tt == "BR": tt = "Br";
      if tt == "MG": tt = "Mg";
      if tt == "CA": tt = "Ca";
      if tt == "MN": tt = "Mn";
      if tt == "FE": tt = "Fe";
      if tt == "ZN": tt = "Zn";
      r = t2r[tt]
      pref = line[:30]
      print('%s%10.4f%10.4f%10.4f %7.4f %6.4f' % (pref, x, y, z, q, r))


parameters = '''
H      2.736     0.0157
C      3.816     0.0860
N      3.648     0.1700
O      3.40      0.2100
S      4.00      0.2500
P      4.20      0.200 
F      3.50      0.0610
Cl     3.896     0.2650
Br     4.44      0.3200
I      4.70      0.4000
Mg     1.30      0.875 
Ca     1.98      0.550 
Na     1.868     0.0027
Mn     1.30      0.875 
Fe     1.30      0.010 
Zn     1.48      0.550 
'''

if __name__ == '__main__':
  if len(sys.argv) == 1:
    print 'Usage:', sys.argv[0], 'PDBQEFile'
  if len(sys.argv) == 2:
    run(parameters, sys.argv[1])
