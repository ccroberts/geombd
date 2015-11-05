#!/usr/bin/env python3

import sys

t2r = {}

for line in open(sys.argv[1], 'r'):
  if line.startswith('atom_par'):
    sp = line.split()
    t = sp[1]
    r = float(sp[2]) / 2.
    t2r[t] = r

for line in open(sys.argv[2], 'r'):
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
