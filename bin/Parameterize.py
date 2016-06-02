#!/usr/bin/env python

import sys, subprocess



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


def process_receptor(af, ff, resname, resid, resdata):
  if resname in ff.keys():
    for line in resdata:
      at = line[12:16].strip()
      at = synonym(at)
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      try:
        param = ff[resname][at]
        #print '%s  1.00  0.00   %8.4f %s' % (line[:54], param[0], param[1])
        print('%s%10.4f%10.4f%10.4f %7.4f %6.4f %s' % (line[:30], x, y, z, param[0], af[param[1]]['r'], param[1])) #PQR
      except KeyError:
        print 'REMARK No charge assign for atom type', at
        #print '%s  1.00  0.00   %8.4f %s' % (line[:54], 0, at[0])
        print('%s%10.4f%10.4f%10.4f %7.4f %6.4f %s' % (line[:30], x, y, z, 0, af[at[0]]['r'], at[0])) #PQR
  else:
    print 'REMARK: No assignment for residue', resname, resid, ' - Using OpenBabel to assign Gasteiger charges.'
    tmpfd = open('/tmp/het.pdb', 'w')
    for line in resdata:
      tmpfd.write(line)
    tmpfd.close()
    subprocess.call(["babel", "-ipdb", "/tmp/het.pdb", "-opdb", "/tmp/heth.pdb"])
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
        #TODO... better atom type parsing
        at = line[12:16].strip()[0]#TODO
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        #print '%s  1.00  0.00   %8.4f %s' % (line[:54], Q[i], at[0]) #PDBQE
        print('%s%10.4f%10.4f%10.4f %7.4f %6.4f %s' % (line[:30], x, y, z, Q[i], af[at]['r'], at)) #PQR
        i += 1



def run(paramfn, pdbfn):
  ff = {}
  af = {}

  for line in open(paramfn, 'r'):
    sp = line.split()
    if len(sp) < 5: continue
    if sp[0] == 'frag_par':
      try:
        ff[sp[1]][sp[2]] = [float(sp[3]), sp[4]]
      except KeyError:
        ff[sp[1]] = {}
        ff[sp[1]][sp[2]] = [float(sp[3]), sp[4]]
    if sp[0] == 'atom_par':
      af[sp[1]] = {'r':float(sp[2])/2., 'eps':float(sp[3]), 'vol':float(sp[4]), 'sol':float(sp[5])}

  resn = None
  resi = -6969
  resd = []
  term = False

  for line in open(pdbfn, 'r'):
    if line.startswith('ATOM') or line.startswith('HETATM'):
      rid = int(line[23:26])
      rnm = line[17:20].strip()
      if rid != resi:
        if resi == -6969:
          tnm = 'N%s' % rnm
          if tnm in ff.keys(): rnm = tnm
        else:
          process_receptor(af, ff, resn, resi, resd)
          resd = []
          if term:
            print 'TER'
        if term:
          tnm = 'N%s' % rnm
          if tnm in ff.keys(): rnm = tnm
          term = False
        resi = rid
        resn = rnm
      resd.append(line)
      if line[12:16].strip() == 'OXT':
        resn = 'C%s' % resn
        term = True
    if line.startswith('END'):
      if len(resd) != 0:
        process_receptor(af, ff, resn, resi, resd)
        resi = -6969
        resd = []
        resn = None
        print 'END'


if __name__ == '__main__':
  if len(sys.argv) == 1:
    print 'Usage:', sys.argv[0], 'ParameterFile', 'PDBFile'
  if len(sys.argv) == 3:
    run(sys.argv[1], sys.argv[2])
