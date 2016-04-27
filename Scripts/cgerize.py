#!/usr/bin/env python

import sys, subprocess




def process_receptor(ff, resname, resid, resdata):
  if resname in ff.keys():
    for line in resdata:
      at = line[12:16].strip()
      try:
        param = ff[resname][at]
        print '%s%s %s  1.00  0.00   %8.4f %s' % (line[:13], param[1], line[15:54], param[0], param[1])
      except KeyError:
        print 'REMARK No charge assign for atom type', at
        print '%s  1.00  0.00   %8.4f %s' % (line[:54], 0, '?')
  else:
    print 'REMARK: No assignment for residue', resname, resid, ' - Using OpenBabel to assign Gasteiger charges.'
    for line in resdata:
      print '%s  1.00  0.00   %8.4f %s' % (line[:54], 0, '?')



def run(params, pdbfn):
  ff = {}

  for line in params.split('\n'):
    sp = line.split()
    if len(sp) < 2: continue
    try:
      ff[sp[0]][sp[1]] = [float(sp[2]), sp[3]]
    except KeyError:
      ff[sp[0]] = {}
      ff[sp[0]][sp[1]] = [float(sp[2]), sp[3]]
    resn = None
    resi = -6969
    resd = []

  for line in open(pdbfn, 'r'):
    if line.startswith('ATOM') or line.startswith('HETATM'):
      rid = int(line[23:26])
      aid = line[13:15]
      rnm = line[17:20].strip()
      if rid != resi:
        if resi == -6969:
          tnm = 'N%s' % rnm
          if tnm in ff.keys(): rnm = tnm
        else:
          process_receptor(ff, resn, resi, resd)
          resd = []
        resi = rid
        resn = rnm
      if aid == 'CA':
        resd.append(line)
      if line[13:16] == 'OXT': resn = 'C%s' % resn
    if line.startswith('END'):
      if len(resd) != 0:
        process_receptor(ff, resn, resi, resd)
        resi = -6969
        resd = []
        resn = None
        print 'END'

parameters = '''
ALA CA  0  a
ARG CA  1  r
ASN CA  0  n
ASP CA -1  d
ASH CA  0  d
CYS CA  0  c
CYM CA -1  c
CYX CA  0  c
GLN CA  0  q
GLU CA -1  e
GLH CA  0  e
GLY CA  0  g
HIS CA  0  h
HID CA  0  h
HIE CA  0  h
HIP CA  0  h
ILE CA  0  i
LEU CA  0  l
LYS CA  1  k
LYN CA  0  k
MET CA  0  m
PHE CA  0  f
PRO CA  0  p
SER CA  0  s
THR CA  0  t
TRP CA  0  w
TYR CA  0  y
VAL CA  0  v
CALA CA  0  a
CARG CA  1  r
CASN CA  0  n
CASP CA -1  d
CASH CA  0  d
CCYS CA  0  c
CCYM CA -1  c
CCYX CA  0  c
CGLN CA  0  q
CGLU CA -1  e
CGLH CA  0  e
CGLY CA  0  g
CHIS CA  0  h
CHID CA  0  h
CHIE CA  0  h
CHIP CA  0  h
CILE CA  0  i
CLEU CA  0  l
CLYS CA  1  k
CLYN CA  0  k
CMET CA  0  m
CPHE CA  0  f
CPRO CA  0  p
CSER CA  0  s
CTHR CA  0  t
CTRP CA  0  w
CTYR CA  0  y
CVAL CA  0  v
NALA CA  0  a
NARG CA  1  r
NASN CA  0  n
NASP CA -1  d
NASH CA  0  d
NCYS CA  0  c
NCYM CA -1  c
NCYX CA  0  c
NGLN CA  0  q
NGLU CA -1  e
NGLH CA  0  e
NGLY CA  0  g
NHIS CA  0  h
NHID CA  0  h
NHIE CA  0  h
NHIP CA  0  h
NILE CA  0  i
NLEU CA  0  l
NLYS CA  1  k
NLYN CA  0  k
NMET CA  0  m
NPHE CA  0  f
NPRO CA  0  p
NSER CA  0  s
NTHR CA  0  t
NTRP CA  0  w
NTYR CA  0  y
NVAL CA  0  v
'''

if __name__ == '__main__':
  if len(sys.argv) == 1:
    print 'Usage:', sys.argv[0], 'PDBFile'
  if len(sys.argv) == 2:
    run(parameters, sys.argv[1])
