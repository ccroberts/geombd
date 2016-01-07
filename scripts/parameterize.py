#!/usr/bin/env python

import sys, subprocess

fffn = sys.argv[1]
pdbfn = sys.argv[2]

ff = {}

for line in open(fffn, 'r'):
  sp = line.split()
  try:
    ff[sp[0]][sp[1]] = [float(sp[2]), sp[3]]
  except KeyError:
    ff[sp[0]] = {}
    ff[sp[0]][sp[1]] = [float(sp[2]), sp[3]]

def synonym(at):
  if at == 'H5\'': return 'H5\'1'
  elif at == 'H5\'\'': return 'H5\'2'
  elif at == 'H2\'': return 'H2\'1'
  elif at == 'H2\'\'': return 'H2\'2'
  elif at == 'HO5\'': return 'HO5\'1'
  elif at == 'HO5\'\'': return 'HO5\'2'
  elif at == 'HO2\'': return 'HO2\'1'
  elif at == 'HO2\'\'': return 'HO2\'2'
  return at


def process_receptor(ff, resname, resid, resdata):
  if resname in ff.keys():
    for line in resdata:
      at = line[12:16].strip()
      at = synonym(at)
      try:
        param = ff[resname][at]
        print '%s  1.00  0.00   %7.3f %s' % (line[:54], param[0], param[1])
      except KeyError:
        print 'REMARK No charge assign for atom type', at
        print '%s  1.00  0.00   %7.3f %s' % (line[:54], 0, at[0])
  else:
    print 'REMARK: No assignment for residue', resname, resid, ' - Using OpenBabel to assign Gasteiger charges.'
    tmpfd = open('/tmp/het.pdb', 'w')
    for line in resdata:
      tmpfd.write(line)
    tmpfd.close()
    subprocess.call(["babel", "-ipdb", "/tmp/het.pdb", "-opdb", "/tmp/heth.pdb", "-h"])
    subprocess.call(["babel", "-ipdb", "/tmp/heth.pdb", "-omol2", "/tmp/het.mol2"])
    Q = []
    tmpfd = open('/tmp/het.mol2', 'r')
    tmpfd.readline()
    tmpfd.readline()
    n = int(tmpfd.readline().split()[0])
    line = tmpfd.readline()
    while not line.startswith("@<TRIPOS>ATOM"):
      line = tmpfd.readline()
    for line in tmpfd:
      if line[0] == '@' or line.strip() == '': break
      sp = line.split()
      q = float(sp[-1])
      Q.append(q)
    i = 0
    for line in open('/tmp/heth.pdb', 'r'):
      if line.startswith('ATOM') or line.startswith('HETATM'):
        at = line[12:16].strip()
        # todo... better atom type parsing
        print '%s  1.00  0.00   %7.3f %s' % (line[:54], Q[i], at[0])
        i += 1




resn = None
resi = -6969
resd = []

for line in open(pdbfn, 'r'):
  if line.startswith('ATOM') or line.startswith('HETATM'):
    rid = int(line[23:26])
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
    resd.append(line)
    if line[13:16] == 'OXT': resn = 'C%s' % resn
  if line.startswith('END'):
    if len(resd) != 0:
      process_receptor(ff, resn, resi, resd)
      resi = -6969
      resd = []
      resn = None
      print 'END'

